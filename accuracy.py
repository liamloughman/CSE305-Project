import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

algorithm_names = {
    1: "sequential_simulation",
    2: "parallel_step_simulation",
    3: "compute_forces_parallel",
    4: "parallel_combined_simulation",
    5: "barnes_hutt_simulation"
}

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

    for algo_id in range(2, 6):
            for body in trajectories:
                if algo_id in trajectories[body]:
                    sequential_df = trajectories[body][1]
                    df = trajectories[body][algo_id]
                    errors = np.sqrt((sequential_df['x'] - df['x'])**2 + (sequential_df['y'] - df['y'])**2)
                    average_errors[algo_id] += errors

            average_errors[algo_id] /= len(bodies)

    plt.figure(figsize=(10, 8))
    for algo_id, errors in average_errors.items():
            plt.plot(errors, label=algorithm_names[algo_id])

    plt.title(f'Average Positional Error per Time Step Relative to {algorithm_names[1]}')
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
    for body in bodies:
        for algorithm in data:
            algo_id = algorithm['algorithm']
            positions = algorithm[body]
            df = pd.DataFrame(positions)
            trajectories[algo_id] = df

        plt.figure(figsize=(10, 8))
        for algo_id, df in trajectories.items():
                x_values = df['x'].to_numpy()
                y_values = df['y'].to_numpy()
                plt.plot(x_values, y_values, label=algorithm_names[algo_id])
        
        plt.title(f'Trajectories of {random_body} by the First 4 Algorithms')
        plt.xlabel('X Position')
        plt.ylabel('Y Position')
        plt.legend()
        plt.savefig(f'combined_trajectory_plot_{body}.png')
        plt.close()

plot_average_error()
plot_single_trajectory()
