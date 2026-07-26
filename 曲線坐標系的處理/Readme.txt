在General interpolatino LBM，jacobian的求解為透過 6oder finite difference 來做求解
以及透過contravriant曲線追蹤求解(r,s)
然而
在 isoparamteric interpolation-based LBM 是透過等參數條件
求解，每一個計算點的jacobian的coefficient，以及透過newton-rapgson method求解(r,s)

而綜合上述，能夠概括兩者的技術個特點為曲線網格系統+semi-lagrange LBM