import matplotlib.pyplot as plt
import numpy as np

from femnodes import Node
from femelement import BaseElement


def draw_structure(is_deformed: bool, elements: list[BaseElement], nodes: list[Node]):
    scaleF = 10

    if is_deformed == 0:
        UnLineColor = [128/255, 128/255, 128/255]
        UnLineThickness = 4
        aLineStyle = '-'
    else:
        UnLineColor = [255/255, 255/255, 255/255]
        UnLineThickness = 3
        aLineStyle = ':'

    if True:
        plt.figure()
        plt.axis('off')

        # Draw members
        for ielem, element in enumerate(elements):
            x_line = []
            y_line = []

            if is_deformed:
                x_lined = []
                y_lined = []

            for end_node in element.get_nodes():
                x_line.append(end_node.get_coordinates()[0])
                y_line.append(end_node.get_coordinates()[1])

                if is_deformed:
                    x_lined.append(end_node.get_coordinates()[0] + scaleF * end_node.get_displacements()[0])
                    y_lined.append(end_node.get_coordinates()[1] + scaleF * end_node.get_displacements()[1])

            plt.plot(x_line, y_line, linestyle=aLineStyle, linewidth=UnLineThickness, color=UnLineColor)

            if is_deformed:
                plt.plot(x_lined, y_lined, linewidth=4, color=[128/255, 128/255, 128/255])

            if not is_deformed:
                xmid = sum(x_line) / 2
                ymid = sum(y_line) / 2
            else:
                xmid = sum(x_lined) / 2
                ymid = sum(y_lined) / 2

            # Label elements
            plt.text(xmid - 0.1, ymid + 0.1, str(ielem + 1), fontsize=14)

        # Draw Nodes
        for inode, node in enumerate(nodes):
            x = node.get_coordinates()[0]
            y = node.get_coordinates()[1]

            if is_deformed:
                xd = x + scaleF * node.get_displacements()[0]
                yd = y + scaleF * node.get_displacements()[1]

            plt.plot(x, y, 'r.', markersize=30)

            if is_deformed:
                plt.plot(xd, yd, 'r.', markersize=30)

            # Label nodes
            if not is_deformed:
                plt.text(x - 0.05, y + 0.05, str(inode + 1), fontsize=14)
            else:
                plt.text(xd - 0.05, yd + 0.05, str(inode + 1), fontsize=14)

        # if not is_deformed:
        #     plt.savefig('theFig.png')
        # else:
        #     plt.savefig('theFigD.png')

        plt.show(block=False)

    return 1
