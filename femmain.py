from femanalysis import *

from input_File_Thema_1_1 import analysis as analysis1
from input_File_Thema_1_2 import analysis as analysis2
from input_File_Thema_1_3 import analysis as analysis3
from input_File_Thema_2 import analysis as analysis4


if __name__ == "__main__":
    print("Thema 1")
    print("Erotima 1")
    print("-"*45)
    analysis1.analyze()
    analysis1.draw_structure()
    analysis1.print_results()
    print("")

    print("Erotima 2")
    print("-"*45)
    analysis2.analyze()
    analysis2.draw_structure()
    analysis2.print_results()
    print("")

    print("Erotima 3")
    print("-"*45)
    analysis3.analyze()
    analysis3.draw_structure()
    analysis3.print_results()
    print("")

    print("-" * 45)
    print("Thema 2")
    print("-" * 45)
    analysis4.analyze()
    analysis4.draw_structure()
    analysis4.print_results()
    print("")

    print("Analysis is Finished.")
    print("-" * 45)
    x = input("Press enter to exit.")
