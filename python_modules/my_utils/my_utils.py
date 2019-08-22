#!/usr/bin/env python3



#============================================
# This module contains various small
# utilities that I'm getting sick and tired
# of programming again and again by hand
#============================================




#====================================
def yesno(*obj):
#====================================
    """
    Ask a yes/no question.
    Returns True if yes, False if no.

    You can use yesno() to print the question the same way
    you would use print() [without the formatting].
    """

    question = " ".join([str(o) for o in obj])

    while True:

        ans = input(question+' [y/n] ')
        try:
            form_ans = ans[0].lower()
        except IndexError:
            # happens for empty input
            form_ans = None

        if form_ans == 'y':
            return True
        elif form_ans == 'n':
            return False
        else:
            print("Didn't understand your input. Try again with 'y' or 'n'.")


    return
