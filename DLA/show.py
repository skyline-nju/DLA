import os
import sys
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


def add_line(ax,
             x_beg,
             y_beg,
             x_end,
             slope,
             label=None,
             xl=None,
             yl=None,
             fontsize="x-large",
             scale="lin",
             c="#7f7f7f"):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    if scale == "lin":
        slope_new = slope * (xmax - xmin) / (ymax - ymin)
    else:
        slope_new = slope * (np.log10(xmax / xmin) / np.log10(ymax / ymin))
    x = np.linspace(x_beg, x_end, 100)
    y = slope_new * (x - x_beg) + y_beg
    ax.plot(x, y, "-.", transform=ax.transAxes, color=c)
    if label is not None:
        width = ax.bbox.width
        height = ax.bbox.height
        deg = np.arctan(slope_new * height / width) * 180 / np.pi
        dx = x_end - x_beg
        if xl is None:
            xl = x_beg + dx * 0.3
        if yl is None:
            yl = y_beg + dx * 0.6 * slope_new
        ax.text(
            xl,
            yl,
            label,
            transform=ax.transAxes,
            rotation=deg,
            color=c,
            fontsize=fontsize)


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


def show_rect(file, a, b, ax=None, fill=False, nmax=None):
    x, y, theta = read(file)
    if ax is None:
        flag_show = True
        ax = plt.subplot(111)
    else:
        flag_show = False
    if nmax is None:
        size = x.size
    else:
        size = nmax
    if fill is False:
        for i in range(size):
            ax.add_patch(
                patches.Rectangle(
                    (x[i], y[i]), a, b, angle=theta[i], fill=fill))
    else:
        c = plt.cm.viridis(np.linspace(0, 1, size))
        for i in range(size):
            ax.add_patch(
                patches.Rectangle(
                    (x[i], y[i]), a, b, angle=theta[i], color=c[i]))

    ax.axis("equal")
    # ax.axis("off")
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
        mass[i] = mass[i - 1] + hist[i]
    # plt.loglog(rbins, mass, "-o")
    # plt.show()
    # plt.close()
    with open("fd_" + file, "w") as f:
        for i in range(mass.size):
            f.write("%f\t%f\n" % (rbins[i], mass[i]))


def plot_r_vs_N():
    def read_file(file):
        with open(file) as f:
            lines = f.readlines()
            r = np.zeros(len(lines))
            N = np.zeros_like(r)
            for i, line in enumerate(lines):
                s = line.replace("\n", "").split("\t")
                r[i] = float(s[0])
                N[i] = float(s[1])
        return r, N

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    for angle in [0, 2, 4, 6, 8, 10]:
        r, N = read_file("fd_%d.dat" % angle)
        ax.loglog(r, N, "o", label=r"$\theta=%d\degree$" % (angle))
    plt.legend(fontsize="x-large")
    plt.xlim(10)
    plt.ylim(1)
    add_line(ax, 0.5, 0.5, 0.9, 1.71, scale="log", label=r"$slope = 1.71$")
    plt.xlabel(r"$r$", fontsize="xx-large")
    plt.ylabel(r"$N$", fontsize="xx-large")
    plt.tight_layout()
    plt.show()
    plt.close()


if __name__ == "__main__":
    # if platform.system() == "Windows":
    #     os.chdir("data")
    #     file = "500_8_120.dat"
    # else:
    #     file = sys.argv[1]
    # fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 7))
    # show_rect(file, 14, 2, ax, fill=True, nmax=None)
    # plt.tight_layout()
    # if platform.system() == "Windows":
    #     plt.show()
    # else:
    #     plt.savefig(file.replace(".dat", ".png"))
    # plt.close()
    # cal_fracal_dimension(file)
    # plot_r_vs_N()

    os.chdir("data/CCW")
    N = 500
    angle = [2, 4, 6, 8, 10, 12]
    tag = ["a", "b", "c", "d", "e", "f"]
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(11, 8))
    for i, ax in enumerate(axes.flat):
        file = "100000_%d_123.dat" % (angle[i])
        show_rect(file, 14, 2, ax, fill=True, nmax=N)
        ax.set_title(r"${\rm (%s)}\ \theta = %d\degree$" % (tag[i], angle[i]))
    plt.suptitle(r"$N=%d$" % (N), fontsize="xx-large")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()
    # plt.savefig("%d.png" % (N))
    plt.close()
