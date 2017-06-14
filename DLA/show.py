import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def show_disk(file):
    with open(file) as f:
        lines = f.readlines()[2:]
        x = np.zeros(len(lines))
        y = np.zeros_like(x)
        for i, line in enumerate(lines):
            s = line.replace("\n", "").split("\t")
            x[i] = float(s[1])
            y[i] = float(s[2])
    ax = plt.subplot(111)
    if x.size <= 10000:
        if x.size <= 5000:
            aa = True
        else:
            aa = False
        for i in range(x.size):
            ax.add_patch(patches.Circle((x[i], y[i]), 1, aa=aa))
    else:
        for i in range(x.size):
            ax.add_patch(patches.CirclePolygon(
                    (x[i], y[i]), 1, resolution=10, aa=False))
    ax.axis("equal")
    plt.show()
    plt.close()


def show_rect(file, a, b):
    with open(file) as f:
        lines = f.readlines()[2:]
        x = np.zeros(len(lines))
        y = np.zeros_like(x)
        theta = np.zeros_like(x)
        for i, line in enumerate(lines):
            s = line.replace("\n", "").split("\t")
            x[i] = float(s[1])
            y[i] = float(s[2])
            theta[i] = float(s[3])
    ax = plt.subplot(111)
    for i in range(x.size):
        ax.add_patch(patches.Rectangle((x[i], y[i]), b, a, angle=theta[i]))
    ax.axis("equal")
    plt.show()
    plt.close()


if __name__ == "__main__":
    ax = plt.subplot(111)
    ax.add_patch(patches.Rectangle((0, 0), 1, 0.5, angle=45))
    ax.add_patch(patches.Rectangle((0, 0), 3.5, 0.5, fill=False))
    ax.add_patch(patches.Rectangle((2, 2), 1, 0.5))
    ax.add_patch(patches.Circle((1, 1), 1))
    ax.axis("equal")
    plt.show()
    # show_disk("N1000.xyz")
