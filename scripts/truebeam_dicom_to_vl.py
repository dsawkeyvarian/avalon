import numpy as np
import pydicom
import argparse
import json
from typing import Union


#Point of this is to make the output more pretty: short objects and mlc_leaf_positions on a single line
class CompactJSONEncoder(json.JSONEncoder):
    """A JSON Encoder that puts small containers on single lines."""

    CONTAINER_TYPES = (list, tuple, dict)
    """Container datatypes include primitives or other containers."""

    MAX_WIDTH = 70
    """Maximum width of a container that might be put on a single line."""

    MAX_ITEMS = 4
    """Maximum number of items in container that might be put on single line."""

    INDENTATION_CHAR = " "

    def __init__(self, *args, **kwargs):
        # using this class without indentation is pointless
        if kwargs.get("indent") is None:
            kwargs.update({"indent": 4})
        super().__init__(*args, **kwargs)
        self.indentation_level = 0

    def encode(self, o):
        """Encode JSON object *o* with respect to single line lists."""
        if isinstance(o, (list, tuple)):
            if self._put_on_single_line(o):
                return "[" + ", ".join(self.encode(el) for el in o) + "]"
            else:
                self.indentation_level += 1
                output = [self.indent_str + self.encode(el) for el in o]
                self.indentation_level -= 1
                return "[\n" + ",\n".join(output) + "\n" + self.indent_str + "]"
        elif isinstance(o, dict):
            if o:
                if self._put_on_single_line(o):
                    return "{ " + ", ".join(f"{self.encode(k)}: {self.encode(el)}" for k, el in o.items()) + " }"
                else:
                    self.indentation_level += 1
                    output = [self.indent_str + f"{json.dumps(k)}: {self.encode(v)}" for k, v in o.items()]
                    self.indentation_level -= 1
                    return "{\n" + ",\n".join(output) + "\n" + self.indent_str + "}"
            else:
                return "{}"
        elif isinstance(o, float):  # Use scientific notation for floats, where appropiate
            return format(o, "g")
        elif isinstance(o, str):  # escape newlines
            o = o.replace("\n", "\\n")
            return f'"{o}"'
        else:
            return json.dumps(o)

    def iterencode(self, o, **kwargs):
        """Required to also work with `json.dump`."""
        return self.encode(o)

    def _put_on_single_line(self, o):
        is_float_list = False
        if isinstance(o, (list, tuple)):
            all_values_floats = True
            for val in o:
                if not isinstance(val, float):
                    all_values_floats = False
            is_float_list = all_values_floats
        
        return is_float_list or (self._primitives_only(o) and len(o) <= self.MAX_ITEMS and len(str(o)) - 2 <= self.MAX_WIDTH)

    def _primitives_only(self, o: Union[list, tuple, dict]):
        if isinstance(o, (list, tuple)):
            return not any(isinstance(el, self.CONTAINER_TYPES) for el in o)
        elif isinstance(o, dict):
            return not any(isinstance(el, self.CONTAINER_TYPES) for el in o.values())

    @property
    def indent_str(self) -> str:
        return self.INDENTATION_CHAR*(self.indentation_level*self.indent)

if __name__ == "__main__":
    import sys
    parser = argparse.ArgumentParser(description="Parses truebeam control points from dicom beam to vl format")
    parser.add_argument('-p', '--plan', help='Dicom Plan to be modified.', required=True)
    parser.add_argument('-b', '--beam', help='Select one beam.', default=0)
    parser.add_argument('-o', '--out', help='Output filename')
    parser.add_argument('-ls', help='List all beams in the plan and their indices', action='store_true')
    args = parser.parse_args()
    plan = pydicom.read_file(args.plan)
    if args.ls == True:
        print(f"Number of beams: {len(plan.BeamSequence)}")
        for idx, beam in enumerate(plan.BeamSequence):
            print(f"Idx: {idx:>2}, Name: {beam.BeamName}")
        sys.exit()
    beam_idx = int(args.beam)

    cumulative_meterset_weights = []
    collimator_angles = []
    gantry_angles = []
    jaw_xs = []
    jaw_ys = []
    leaf_sequences = []
    num_cps = 0
    print(f"Selected Beam idx : {beam_idx}")
    print(f"Selected Beam Name: { plan.BeamSequence[beam_idx].BeamName}")
    for cp in plan.BeamSequence[beam_idx].ControlPointSequence:
        num_cps += 1
        cumulative_meterset_weights.append(float(cp.CumulativeMetersetWeight))
        gantry_angles.append(float(cp.GantryAngle) if "GantryAngle" in cp else gantry_angles[-1])
        collimator_angles.append(float(cp.BeamLimitingDeviceAngle) if "BeamLimitingDeviceAngle" in cp else collimator_angles[-1])
        blds = {}
        if "BeamLimitingDevicePositionSequence" in cp:
            for ld_values in cp.BeamLimitingDevicePositionSequence:
                if ld_values.RTBeamLimitingDeviceType=='MLCX':
                    blds['MLCX'] = np.array(ld_values.LeafJawPositions)
                elif ld_values.RTBeamLimitingDeviceType=='X' or ld_values.RTBeamLimitingDeviceType=='ASYMX':
                    jaw_x_pos = np.array(ld_values.LeafJawPositions)
                    assert len(jaw_x_pos) == 2 
                    blds['X'] = jaw_x_pos
                elif ld_values.RTBeamLimitingDeviceType=='Y' or ld_values.RTBeamLimitingDeviceType=='ASYMY':
                    jaw_y_pos = np.array(ld_values.LeafJawPositions)
                    assert len(jaw_y_pos) == 2 
                    blds['Y'] = jaw_y_pos
                else:
                    raise RuntimeError(f"Unexpected beam limiting device type {ld_values.RTBeamLimitingDeviceType}")
        leaf_sequences.append(blds['MLCX'] if 'MLCX' in blds else np.copy(leaf_sequences[-1]))
        jaw_xs.append(blds['X'] if 'X' in blds else np.copy(jaw_xs[-1]))
        jaw_ys.append(blds['Y'] if 'Y' in blds else np.copy(jaw_ys[-1]))

    assert len(collimator_angles) == num_cps
    assert len(gantry_angles) == num_cps
    assert len(jaw_xs) == num_cps
    assert len(jaw_ys) == num_cps

    #Assemble .json plan
    num_leaves_in_bank = 60 #Hardcoded truebeam hd/millennium120 value
    plan = []
    for i in range(num_cps):
        control_point = {
            "cumulative_meterset_weight": cumulative_meterset_weights[i],
            "angles": {"unit": "deg", "gantry": gantry_angles[i], "collimator": collimator_angles[i]},
            "jaw_positions": {
                "unit": "mm",
                "x1": jaw_xs[i][0], "x2": jaw_xs[i][1],
                "y1": jaw_ys[i][0], "y2": jaw_ys[i][1]
            },
            "mlc_leaf_positions": {
                "unit": "mm",
                "bank_X1": [val for val in leaf_sequences[i][0:num_leaves_in_bank]],
                "bank_X2": [val for val in leaf_sequences[i][num_leaves_in_bank:]]
            }
        }
        plan.append(control_point)
    if not args.out:
        print("No output filename given, output not written")
    else:
        with open(args.out, "w", encoding="utf-8") as outf:
            json.dump({"control_points": plan}, outf, indent=2, cls=CompactJSONEncoder)
