import pandas as pd
import seaborn as sns
import  matplotlib.pyplot as plt

def main():
    func = pd.read_csv("func_inspect.csv")
    print(func.head())
    plt.plot(func.iloc[:,0], func.iloc[:,1])
    plt.show()

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
