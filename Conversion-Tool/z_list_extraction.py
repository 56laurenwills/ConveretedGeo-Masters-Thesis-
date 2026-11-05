import re
from collections import defaultdict

# Name_map = {
#     (1, 1): "neg_disk_envelope",
#     (1, 2): "pos_disk_envelope",
#     (4, 1): "innerPixel_neg-ec_right",
#     (4, 2): "innerPixel_neg-ec_left",
#     (5, 1): "innerPixel_barrel",
#     (5, 2): "innerpixel_other-barrel",
#     (6, 1): "innerpixel_pos-ec_right",
#     (6, 2): "innerpixel_neg-ec_right",
#     (9, 1): "outerPixel_neg_ec_right",
#     (9, 2): "outerPixel_neg_ec_left",
#     (10, 1): "outerPixel_barrel_1",
#     (10, 2): "outerPixel_barrel_2",
#     (19, 1): "strip_ec_right",
#     (19, 2): "strip_ec_left",
#     (20, 1): "neg_strip_envelope", 
#     (20, 2): "pos_strip_envelope",

# }


def generate_config(filepath, geometry_tag): #name_map
    direct_lookup = {}

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tokens = [t.strip() for t in line.split("|")]
            #if len(tokens) != 5:
            #    continue
            if len(tokens) == 5:
                vol_tok, bnd_tok, lay_tok, volume_name, z_str = tokens
                vol = int(vol_tok.split("=")[1])
                bnd = int(bnd_tok.split("=")[1])
                lay = int(lay_tok.split("=")[1])
                z = abs(float(z_str))

            if len(tokens) == 4:
                vol_tok, bnd_tok, volume_name, z_str = tokens
                vol = int(vol_tok.split("=")[1])
                bnd = int(bnd_tok.split("=")[1])
                lay = 0
                z = abs(float(z_str))

            if vol not in (1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 19, 20, 21): #if (vol, bnd) not in name_map and
                continue

            direct_lookup[(vol, bnd, lay)] = z

    print(f'"{geometry_tag}": {{')
    print('    "z_lookup": {')
    for (vol, bnd, lay), z in sorted(direct_lookup.items()):
        print(f'        ({vol}, {bnd}, {lay}): {z:.3f},')
    print("    },")
    print("},")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", default="disc_portal_z_values.txt")
    parser.add_argument("--tag", default="Run4-04-00")
    args = parser.parse_args()
    generate_config(args.file, args.tag) #name_map=Name_map)

