import pandas as pd
import seaborn as sns
import  matplotlib.pyplot as plt

def main():
    df = pd.read_csv("Gauss_Legendre_table.txt", header = 1)
    df = df.rename(columns = lambda x: x.strip())
    print(df.head())
    print(df.dtypes)
    df[["Result with Gauss-Legendre", "Exact result", "Relative error"]].plot()
    plt.show()

    df["Calculation time [s]"].plot()
    plt.show()
    
if __name__ == "__main__":
    main()
