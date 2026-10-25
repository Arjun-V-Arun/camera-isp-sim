"""Convert any photo (JPG, PNG, ...) to an 8-bit PPM that `isp_sim make-raw` can use as ground truth.

Usage:
    python tools/image_to_ppm.py photo.jpg photo.ppm [--max-side 2048]

Needs Pillow (pip install pillow). Dimensions are trimmed to even numbers so the Bayer tile fits exactly.
"""

import argparse

from PIL import Image


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("src")
    parser.add_argument("dst")
    parser.add_argument("--max-side", type=int, default=0, help="downscale so the longer side is at most this")
    args = parser.parse_args()

    img = Image.open(args.src).convert("RGB")
    if args.max_side and max(img.size) > args.max_side:
        scale = args.max_side / max(img.size)
        img = img.resize((round(img.width * scale), round(img.height * scale)), Image.LANCZOS)
    img = img.crop((0, 0, img.width - img.width % 2, img.height - img.height % 2))
    img.save(args.dst, format="PPM")
    print(f"Wrote {args.dst} ({img.width}x{img.height})")


if __name__ == "__main__":
    main()
