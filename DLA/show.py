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
            ax.add_patch(
                patches.CirclePolygon(
                    (x[i], y[i]), 1, resolution=10, aa=False))
    ax.axis("equal")
    plt.show()
    plt.close()


def show_rect(file, a, b, ax=None, fill=False):
    with open(file) as f:
        lines = f.readlines()
        x = np.zeros(len(lines))
        y = np.zeros_like(x)
        theta = np.zeros_like(x)
        for i, line in enumerate(lines):
            s = line.replace("\n", "").split("\t")
            x[i] = float(s[0])
            y[i] = float(s[1])
            theta[i] = float(s[2])
    if ax is None:
        flag_show = True
        ax = plt.subplot(111)
    else:
        flag_show = False
    for i in range(x.size):
        ax.add_patch(
            patches.Rectangle((x[i], y[i]), a, b, angle=theta[i], fill=fill))
    ax.axis("equal")
    if flag_show:
        plt.show()
        plt.close()


if __name__ == "__main__":
    ax = plt.subplot(111)
    show_rect("rect_50.dat", 14, 2, ax)
    # show_rect("traj.dat", 14, 2, ax, fill=False)
    plt.show()
    plt.close()
    # with open("traj.dat") as f:
    #     lines = f.readlines()
    #     x = np.zeros(len(lines))
    #     y = np.zeros_like(x)
    #     for i, line in enumerate(lines):
    #         s = line.split("\t")
    #         x[i] = float(s[0])
    #         y[i] = float(s[1])
    # plt.plot(x, y, ".")
    # plt.show()
