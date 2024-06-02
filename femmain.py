from femanalysis import *

from input_File_Thema_1_1 import analysis as analysis1
from input_File_Thema_1_2 import analysis as analysis2
from input_File_Thema_1_3 import analysis as analysis3
from input_File_Thema_2 import analysis as analysis4


if __name__ == "__main__":
    print("Erotima 1")
    print("-"*45)
    analysis1.analyze(True)
    analysis1.print_results()
    print("")

    print("Erotima 2")
    print("-"*45)
    analysis2.analyze(True)
    analysis2.print_results()
    print("")

    print("Erotima 3")
    print("-"*45)
    analysis3.analyze(True)
    analysis3.print_results()
    print("")

    print("Thema 2")
    print("-" * 45)
    analysis4.analyze(True)
    print("")

    print("Analysis is Finished.")
    print("-"*45)
    x = input("Press enter to exit.")
