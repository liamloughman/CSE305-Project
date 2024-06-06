import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

def plot_average_error():
    with open('trajectory_results.json', 'r') as file:
        data = json.load(file)

    first_algo_data = data[0]
    bodies = [key for key in first_algo_data.keys() if key.startswith('body')]

    trajectories = {body: {} for body in bodies}

    for algorithm in data:
        algo_id = algorithm['algorithm']
        for body in bodies:
            positions = algorithm[body]
            df = pd.DataFrame(positions)
            trajectories[body][algo_id] = df

    steps_count = len(trajectories[bodies[0]][1])
    average_errors = {algo_id: np.zeros(steps_count) for algo_id in range(2, 6)}

    for algo_id in range(2, 4):
            for body in trajectories:
                if algo_id in trajectories[body]:
                    sequential_df = trajectories[body][1]
                    df = trajectories[body][algo_id]
                    errors = np.sqrt((sequential_df['x'] - df['x'])**2 + (sequential_df['y'] - df['y'])**2)
                    average_errors[algo_id] += errors

            average_errors[algo_id] /= len(bodies)

    plt.figure(figsize=(10, 8))
    for algo_id, errors in average_errors.items():
        if algo_id != 5:
            plt.plot(errors, label=f'Algorithm {algo_id}')

    plt.title('Average Positional Error per Time Step Relative to Algorithm 1')
    plt.xlabel('Time Step')
    plt.ylabel('Average Positional Error')
    plt.legend()
    plt.savefig('time_step_average_error_comparison.png')
    plt.close()

def plot_single_trajectory():
    with open('trajectory_results.json', 'r') as file:
        data = json.load(file)

    bodies = [key for key in data[0].keys() if key.startswith('body')]
    random_body = random.choice(bodies)

    trajectories = {}
    for algorithm in data:
        algo_id = algorithm['algorithm']
        positions = algorithm[random_body]
        df = pd.DataFrame(positions)
        trajectories[algo_id] = df

    plt.figure(figsize=(10, 8))
    for algo_id, df in trajectories.items():
        x_values = df['x'].to_numpy()
        y_values = df['y'].to_numpy()
        plt.plot(x_values, y_values, label=f'Algorithm {algo_id}')
    plt.title(f'Trajectories of {random_body} by All Algorithms')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.legend()
    plt.savefig('combined_trajectory_plot.png')
    plt.close()

    sequential_df = trajectories[1]
    
    plt.figure(figsize=(10, 8))
    for algo_id, df in trajectories.items():
        if algo_id != 1:
            error = np.sqrt((sequential_df['x'].to_numpy() - df['x'].to_numpy())**2 + (sequential_df['y'].to_numpy() - df['y'].to_numpy())**2)
            plt.plot(error, label=f'Error between Algo 1 and Algo {algo_id}')

    plt.title(f'Error Between Sequential (Algo 1) and Other Algorithms')
    plt.xlabel('Time Step')
    plt.ylabel('Positional Error')
    plt.legend()
    plt.savefig('error_comparison_plot.png')
    plt.close()

plot_average_error()
plot_single_trajectory()