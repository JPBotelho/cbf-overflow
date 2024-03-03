import pandas as pd
import matplotlib.pyplot as plt

available_styles = plt.style.available
print(available_styles)
plt.style.use("seaborn-v0_8-paper")
# Load the CSV file into a DataFrame
df = pd.read_csv('data.csv')

# Extract the desired columns
columns_to_extract = ['eX', 'tU']

data = df[columns_to_extract]

# Plot the data
plt.figure(figsize=(10, 6))
markers = ["D", "8", "D"]

labels = ["Exact Value", "Upper bound 2", "Upper Bound 1"]
lineStyles = ["solid", "dashed", "dotted"]
for i, column in enumerate(data.columns):
    plt.plot(data.index[2408:2414], data[column].iloc[2408:2414], label=labels[i],
             ls=lineStyles[i], marker=markers[i], markersize=10)

# plt.xlabel('n')
# plt.ylabel('Overflow rate')
plt.yscale('log')
# plt.title('Overflow Rate of Counting Bloom Filters')
# plt.legend()
plt.grid(False)
plt.show()
