import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("./data/results.csv")

# log time
df["log_time"] = np.log10(df["Average Time (s)"])

# Find Pareto front
def pareto_front(data):
    points = data[["log_time", "Sum MAE"]].values

    pareto = []

    for i, point in enumerate(points):
        dominated = False
        
        for j, other in enumerate(points):
            if i == j:
                continue

            # other is better in both
            if (
                other[0] <= point[0]
                and other[1] <= point[1]
                and (
                    other[0] < point[0]
                    or other[1] < point[1]
                )
            ):
                dominated = True
                break

        if not dominated:
            pareto.append(i)

    return pareto


front = pareto_front(df)

plt.figure(figsize=(8,6))


pareto_df = df.iloc[front].sort_values("log_time")

label_names = {
    "MLIP": "MLIP",
    "SE": "Semi-empirical",
    "DFT": "DFT"
}

colors = {
    "MLIP": "blue",
    "SE": "green",
    "DFT": "purple"
}

for category, group in df.groupby("Type"):

    plt.scatter(
        group["log_time"],
        group["Sum MAE"],
        s=60,
        alpha=0.7,
        label=label_names.get(category, category),
        color=colors.get(category, "gray")
    )

# Connect Pareto points
plt.plot(
    pareto_df["log_time"],
    pareto_df["Sum MAE"],
    color="red",
    linestyle="--",
    linewidth=1.5,
    zorder=2
)

# Specific methods to label
labels = [
    "eSEN-MD-direct-all-OMOL",
    "Orb-v3-direct-inf-OMAT",
    "ωB97M-V",
    "Orbmol-v1-direct"
]

for method in labels:
    row = df[df["Method"] == method]

    if not row.empty:
        plt.annotate(
            method,
            (row["log_time"].values[0], row["Sum MAE"].values[0]),
            xytext=(10, 10),
            textcoords="offset points",
            fontsize=9,
            arrowprops=dict(
                arrowstyle="->",
                linewidth=0.8
            )
        )

from numpy.polynomial import polynomial as poly


plt.xlabel("log10(Average Time (s))")
plt.ylabel("Sum MAE")
plt.title("TRIP50 Accuracy vs Computational Cost")

plt.legend()
plt.grid(True)

plt.savefig(
    "images/pareto_front.png",
    dpi=300,
    bbox_inches="tight"
)