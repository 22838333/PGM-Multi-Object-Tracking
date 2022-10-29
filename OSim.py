# from re import X
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

class object:
    def __init__(self, tag, xpos, ypos):
        self.tag = tag
        self.direction = [0,0]
        self.XPosList = []
        self.YPosList = []
        self.XPosList.append(xpos)
        self.YPosList.append(ypos)

    def move(self, steps, velo):
        if self.direction == [0,0]:
            rad = npr.rand()*2*np.pi
            self.direction = [np.cos(rad), np.sin(rad)]
        for i in range(steps):
            self.XPosList.append(self.XPosList[-1]+self.direction[0]*velo)
            self.YPosList.append(self.YPosList[-1]+self.direction[1]*velo)

    def accelerate(self, steps, velo, accel):
        if self.direction == [0,0]:
            rad = npr.rand()*2*np.pi
            self.direction = [np.cos(rad), np.sin(rad)]
        for i in range(steps):
            velo += accel
            self.XPosList.append(self.XPosList[-1]+self.direction[0]*velo)
            self.YPosList.append(self.YPosList[-1]+self.direction[1]*velo)

    def curve(self, steps, velo, accel):
        if self.direction == [0,0]:
            rad = npr.rand()*2*np.pi
            self.direction = [np.cos(rad)*velo, np.sin(rad)*velo]
            rad = npr.rand()*2*np.pi
            self.acc = [np.cos(rad)*accel, np.sin(rad)*accel]
        for i in range(steps):
            self.direction = [self.direction[0]+self.acc[0], self.direction[1]+self.acc[1]]
            self.XPosList.append(self.XPosList[-1]+self.direction[0])
            self.YPosList.append(self.YPosList[-1]+self.direction[1])

    def measure(self):
        measXPosList = self.XPosList + npr.normal(size = len(self.XPosList))
        measYPosList = self.YPosList + npr.normal(size = len(self.YPosList))
        return [measXPosList,measYPosList]

    def drivel(self, size):
        xStuff = np.zeros_like(self.XPosList) + size*npr.normal(size = len(self.XPosList))
        yStuff = np.zeros_like(self.YPosList) + size*npr.normal(size = len(self.YPosList))
        return [xStuff,yStuff]

    def abnormalDrivel(self, size):
        xStuff = np.zeros_like(self.XPosList) + size*2*(npr.random(size = len(self.XPosList))-0.5)
        yStuff = np.zeros_like(self.YPosList) + size*2*(npr.random(size = len(self.YPosList))-0.5)
        return [xStuff,yStuff]

    def dirAbnormalDrivel(self, size):
        xStuff = self.XPosList + size*2*(npr.random(size = len(self.XPosList))-0.5)
        yStuff = self.YPosList + size*2*(npr.random(size = len(self.YPosList))-0.5)
        return [xStuff,yStuff]

    def clutterMeasure(self, clutterProb, clutterOffset):
        xclutt = np.zeros_like(self.XPosList)
        yclutt = np.zeros_like(self.YPosList)
        for i in range(len(self.XPosList)):
            if(npr.rand()<clutterProb and i != 0):
                # xclutt[i] = (npr.rand()-0.5)*2*0.5*clutterOffset
                # xclutt[i] += np.sign(xclutt[i])*0.5*clutterOffset
                # yclutt[i] = (npr.rand()-0.5)*2*0.5*clutterOffset
                # yclutt[i] += np.sign(yclutt[i])*0.5*clutterOffset
                xclutt[i] = (npr.rand()-0.5)*2*clutterOffset
                yclutt[i] = (npr.rand()-0.5)*2*clutterOffset

        measXPosList = self.XPosList + npr.normal(size = len(self.XPosList)) + xclutt
        measYPosList = self.YPosList + npr.normal(size = len(self.YPosList)) + yclutt


        return [measXPosList,measYPosList]

    # def clutterMeasure(self, clutterProb, clutterOffset):
    #     xclutt = np.zeros_like(self.XPosList)
    #     yclutt = np.zeros_like(self.YPosList)
    #     for i in range(len(self.XPosList)):
    #         if(npr.rand()<clutterProb):
    #             xclutt[i] = (npr.rand()-0.5)*0.5*clutterOffset
    #             xclutt[i] += np.sign(xclutt[i])*0.5*clutterOffset
    #         if(npr.rand()<clutterProb):
    #             yclutt[i] = (npr.rand()-0.5)**clutterOffset
    #             yclutt[i] += np.sign(yclutt[i])*0.5*clutterOffset

    #     measXPosList = self.XPosList + npr.normal(size = len(self.XPosList)) + xclutt
    #     measYPosList = self.YPosList + npr.normal(size = len(self.YPosList)) + yclutt


    #     return [measXPosList,measYPosList]


def main():
    o2 = object("2", 16, 10)
    print(o2.tag, " should now exist with x at: ", o2.XPosList[0])
    o2.move(3,2)
    plt.plot(o2.XPosList, o2.YPosList, 'bo')

# main()