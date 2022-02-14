import sys
from step_1_browser import setup_browser

if __name__ == "__main__":
    try:
        setup_browser(True)
    except:
        print("Script fell down, but it doesn't matter")
