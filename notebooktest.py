import Tkinter as tk
import ttk
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends import backend_tkagg



root = tk.Tk()

mygreen = "#74B096"
myred = "#F89B34"

style = ttk.Style()

style.theme_create( "glycosylator", parent="alt", settings={
        "TNotebook": {"configure": {"tabmargins": [2, 5, 2, 0], "background": mygreen } },
        "TNotebook.Tab": {
            "configure": {"padding": [5, 1], "background": mygreen },
            "map":       {"background": [("selected", myred)],
                          "expand": [("selected", [1, 1, 1, 0])] } }, 
        "TFrame": {"configure": {"background": mygreen}},
        "TButton": {"configure": {"background": mygreen}}                  
         } )

style.theme_use("glycosylator")

note = ttk.Notebook(root)
info = ttk.Frame(note, width=300, height=200)
note.add(info, text = 'Info')
common_glycans = ttk.Frame(note, width=300, height=200)
note.add(common_glycans, text = 'Common glycans')
note.pack(expand=1, fill='both', padx=5, pady=5)
my_glycans = ttk.Frame(note, width=300, height=200)
note.add(my_glycans, text = 'My glycans')
note.pack(expand=1, fill='both', padx=5, pady=5)

ttk.Button(common_glycans, text='glycosylate!').pack(fill='x')

root.mainloop()