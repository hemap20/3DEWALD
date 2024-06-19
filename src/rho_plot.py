import matplotlib.pyplot as plt
import pandas as pd

# Read the CSV file
df = pd.read_csv('rho_output.csv')

# Extract data for plotting
rho_values = df['rho']
avg_energy = df['Average Energy per Particle']
avg_pressure = df['Average Pressure']
excess_chem_potential = df['Excess Chemical Potential']

# Create plots
plt.figure(figsize=(12, 8))

# Average Energy per Particle plot
plt.subplot(2, 2, 1)
plt.plot(rho_values, avg_energy, marker='o', linestyle='-', color='b')
plt.xlabel('Density (rho)')
plt.ylabel('Average Energy per Particle')
plt.title('Average Energy per Particle vs Density')

# Average Pressure plot
plt.subplot(2, 2, 2)
plt.plot(rho_values, avg_pressure, marker='o', linestyle='-', color='r')
plt.xlabel('Density (rho)')
plt.ylabel('Average Pressure')
plt.title('Average Pressure vs Density')

# Excess Chemical Potential plot
plt.subplot(2, 2, 3)
plt.plot(rho_values, excess_chem_potential, marker='o', linestyle='-', color='g')
plt.xlabel('Density (rho)')
plt.ylabel('Excess Chemical Potential')
plt.title('Excess Chemical Potential vs Density')

# Show plots
plt.tight_layout()
plt.show()
