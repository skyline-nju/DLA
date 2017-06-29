import os
import numpy as np
import matplotlib
import platform
if platform.system() == "Windows":
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
else:
    matplotlib.use("Agg")
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


def read(file):
    with open(file) as f:
        lines = f.readlines()
        x = np.zeros(len(lines))
        y = np.zeros_like(x)
        theta = np.zeros_like(x)
        for i, line in enumerate(lines):
            s = line.replace("\n", "").split("\t")
            x[i] = float(s[1])
            y[i] = float(s[2])
            theta[i] = float(s[3])
    return x, y, theta


def show_rect(file, a, b, ax=None, fill=False):
    x, y, theta = read(file)
    if ax is None:
        flag_show = True
        ax = plt.subplot(111)
    else:
        flag_show = False
    if fill is False:
        for i in range(x.size):
            ax.add_patch(
                patches.Rectangle(
                    (x[i], y[i]), a, b, angle=theta[i], fill=fill))
    else:
        c = plt.cm.viridis(np.linspace(0, 1, x.size))
        for i in range(x.size):
            ax.add_patch(
                patches.Rectangle(
                    (x[i], y[i]), a, b, angle=theta[i], color=c[i]))

    ax.axis("equal")
    ax.axis("off")
    if flag_show:
        plt.show()
        plt.close()


def cal_fracal_dimension(file):
    x, y, theta = read(file)
    r = np.sqrt(x * x + y * y)
    rmax = r.max()
    print(rmax)
    rbins = np.logspace(0, np.log2(rmax), 50, base=2)
    bins = np.zeros(rbins.size + 1)
    bins[0] = 0
    bins[1:] = rbins
    hist, bin_edge = np.histogram(r, bins=bins)
    mass = np.zeros_like(rbins)
    mass[0] = hist[0]
    for i in range(1, mass.size):
        mass[i] = mass[i-1] + hist[i]
    plt.loglog(rbins, mass, "-o")
    plt.show()
    plt.close()
    for i in range(mass.size):
        print(rbins[i], mass[i])


if __name__ == "__main__":
    if platform.system() is "Windows":
        os.chdir("data")
    file = "100000_2_2.dat"
    # fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    # show_rect(file, 14, 2, ax, fill=True)

    # plt.tight_layout()
    # if sys.platform == "Windows":
    #     plt.show()
    # else:
    #     plt.savefig(file.replace(".dat", ".png"))
    # plt.close()
    cal_fracal_dimension(file)
