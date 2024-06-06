import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_all_results(directory):
    algorithm_names = {
        1: "sequential_simulation",
        2: "parallel_step_simulation",
        3: "parallel_distinc_simulation",
        4: "parallel_combined_simulation",
        5: "barnes_hutt_simulation"
    }
    
    for n in range(1, 1001):
        file_path = os.path.join(directory, f"results{n}.json")
        if os.path.exists(file_path):
            df = pd.read_json(file_path)
            plt.figure(figsize=(10, 6))
            for algorithm in sorted(df['algorithm'].unique()):
                subset = df[df['algorithm'] == algorithm]
                threads = subset['threads'].values
                times = subset['time'].values
                plt.plot(threads, times, label=f'{algorithm_names[algorithm]}')
            plt.xlabel('Number of Threads')
            plt.ylabel('Time (seconds)')
            plt.title(f'Performance of Algorithms for {n} bodies.')
            plt.legend()
            plt.grid(True)
            plt.savefig(os.path.join(directory, f'plot_results{n}.png'))

plot_all_results('')
