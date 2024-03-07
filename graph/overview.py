import pandas as pd
import matplotlib.pyplot as plt

available_styles = plt.style.available
print(available_styles)
plt.style.use("seaborn-v0_8-paper")
# Load the CSV file into a DataFrame
df = pd.read_csv('data.csv')

diff = df['tU'] - df['eX']

# Extract the desired columns
columns_to_extract = ['eX', 'tU', 'nU']

data = df[columns_to_extract]

# Plot the data
plt.figure(figsize=(10, 6))

labels = ["Exact Value", "Upper Bound (2)", "Upper Bound (1)"]
markers = ["D", "8", "D"]
lineStyles = ["solid", "solid", "solid"]
lineWidth = [2, 7, 2]
i = 0
for column in data.columns:
    plt.plot(data.index, data[column], label=labels[i],
             ls=lineStyles[i], linewidth=lineWidth[i], zorder=3-i)  # , marker=markers[i])
    i += 1

plt.grid(False)

plt.xlabel('n')
plt.ylabel('Overflow rate')
plt.yscale('log')
plt.legend()
plt.show()
