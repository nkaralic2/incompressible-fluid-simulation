"""
visualize.py — Visualization for the PIC/FLIP fluid simulation.

Reads particle CSV files from the output directory and creates
either a static plot, an animation, or exports an MP4/GIF.

Usage:
    python visualize.py                          # Animate from ./output/
    python visualize.py --dir output --save animation.mp4
    python visualize.py --frame 100              # Show single frame
    python visualize.py --dir output --save animation.gif --fps 30
"""

import os
import sys
import glob
import re
import argparse
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend for headless use
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
except ImportError:
    print("Error: matplotlib is required. Install with: pip install matplotlib")
    sys.exit(1)


def load_particles(filepath):
    """Load particle data from a CSV file."""
    try:
        if os.path.getsize(filepath) == 0:
            return None
        data = np.loadtxt(filepath, delimiter=',', skiprows=1)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None


def find_csv_files(directory):
    """Find and sort all particle CSV files in a directory."""
    pattern = os.path.join(directory, "particles_*.csv")
    files = sorted(glob.glob(pattern), key=frame_number)
    return files


def frame_number(filepath):
    """Extract the numeric simulation step from a particle filename."""
    match = re.search(r"particles_(\d+)\.csv$", os.path.basename(filepath))
    return int(match.group(1)) if match else -1

def plot_frame(data, step, domain_width=4.0, domain_height=4.0, ax=None):
    """Plot a single frame of particle positions."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    ax.clear()
    if data is not None and len(data) > 0:
        x = data[:, 0]
        y = data[:, 1]
        # Color by velocity magnitude for visual interest
        if data.shape[1] >= 4:
            vel = np.sqrt(data[:, 2]**2 + data[:, 3]**2)
            scatter = ax.scatter(x, y, s=0.5, c=vel, cmap='coolwarm',
                               vmin=0, vmax=5, alpha=0.8)
        else:
            ax.scatter(x, y, s=0.5, c='blue', alpha=0.8)

    ax.set_xlim(0, domain_width)
    ax.set_ylim(0, domain_height)
    ax.set_aspect('equal')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title(f'FLIP Fluid Simulation — Step {step}')
    ax.set_facecolor('#1a1a2e')

    return ax


def create_animation(files, domain_width=4.0, domain_height=4.0, fps=30):
    """Create a matplotlib animation from CSV files."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    def update(frame_idx):
        filepath = files[frame_idx]
        step = frame_number(filepath)
        data = load_particles(filepath)
        plot_frame(data, step, domain_width, domain_height, ax)
        return []

    anim = animation.FuncAnimation(fig, update, frames=len(files),
                                    interval=1000//fps, blit=False)
    return fig, anim


def main():
    parser = argparse.ArgumentParser(description='FLIP Fluid Simulation Visualizer')
    parser.add_argument('--dir', default='output', help='Directory with particle CSV files')
    parser.add_argument('--frame', type=int, default=None, help='Show a single frame')
    parser.add_argument('--save', default=None, help='Save animation to file (mp4/gif/png)')
    parser.add_argument('--fps', type=int, default=30, help='Frames per second for animation')
    parser.add_argument('--width', type=float, default=4.0, help='Domain width')
    parser.add_argument('--height', type=float, default=4.0, help='Domain height')
    args = parser.parse_args()

    files = find_csv_files(args.dir)
    if not files:
        print(f"No particle CSV files found in '{args.dir}'")
        sys.exit(1)

    print(f"Found {len(files)} frames in '{args.dir}'")

    if args.frame is not None:
        # Show single frame
        target = None
        for f in files:
            step = frame_number(f)
            if step == args.frame:
                target = f
                break
        if target is None:
            print(f"Frame {args.frame} not found")
            sys.exit(1)
        data = load_particles(target)
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        plot_frame(data, args.frame, args.width, args.height, ax)
        if args.save:
            fig.savefig(args.save, dpi=150, bbox_inches='tight')
            print(f"Saved to {args.save}")
        else:
            plt.show()
    else:
        # Create animation
        fig, anim = create_animation(files, args.width, args.height, args.fps)

        if args.save:
            ext = os.path.splitext(args.save)[1].lower()
            if ext == '.gif':
                anim.save(args.save, writer='pillow', fps=args.fps)
            elif ext == '.mp4':
                anim.save(args.save, writer='ffmpeg', fps=args.fps)
            else:
                anim.save(args.save, fps=args.fps)
            print(f"Animation saved to {args.save}")
        else:
            plt.show()


if __name__ == '__main__':
    main()
