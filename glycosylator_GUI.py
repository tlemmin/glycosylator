import Tkinter as tk
import ttk
class GlycosylatorGUI(tk.Tk):

   def __init__(self):
       tk.Tk.__init__(self)
       self.title('Glycosylator')
       self.geometry('500x700')
       #self.resizable(False, False)
       
       # Create menubar
       self.menubar = tk.Menu(master=self, bg="lightgrey", fg="black")
       # file menu
       self.file_menu = tk.Menu(self.menubar, tearoff=0, bg="lightgrey", fg="black") 
       self.file_menu.add_command(label="Open glycoprotein", accelerator = "Ctrl+O") 
       self.file_menu.add_command(label="Save glycoprotein", accelerator = "Ctrl+S") 
       self.menubar.add_cascade(label="File", menu=self.file_menu)
       self.config(menu=self.menubar)
       # glycan menu
       self.gl_menu = tk.Menu(self.menubar, tearoff=0, bg="lightgrey", fg="black") 
       self.gl_menu.add_command(label="Import glycan library", accelerator = "Ctrl+I") 
       self.gl_menu.add_command(label="Export glycan library", accelerator = "Ctrl+E")
       self.menubar.add_cascade(label="Glycan library", menu=self.gl_menu)

       self.config(menu=self.menubar)
       # Create and layout main frames
       self.right_frame = tk.Frame(self, width= 250, height= 500, bg='red')
       self.left_frame = tk.Frame(self, width = 250, height = 500, bg = 'green')
       self.left_frame.grid(column=0, row=0)
       self.right_frame.grid(column=1, row=0)
       
       # Create widget for right frame       
       self.glycosylator_logo = tk.PhotoImage(file="glycosylator_logo.gif")
       self.w_logo = tk.Label(self.right_frame, image=self.glycosylator_logo)
       self.chain_label = tk.Label(self.right_frame, text="Chain")
       options = []
       self.chain = tk.StringVar(self.right_frame)
#       self.chain_menu = tk.OptionMenu(self.right_frame, self.chain, *options)
       self.sequon_label = tk.Label(self.right_frame, text="Sequon")
       self.sequon = tk.StringVar(self.right_frame)
#       self.sequon_menu = tk.OptionMenu(self.right_frame, self.sequon, *options)
       self.glycan_label = tk.Label(self.right_frame, text="Glycan")
       self.glycan_2D = tk.Frame(self.right_frame, width = 150, height = 150, bg = 'cyan') 
       self.glycosylate = tk.Button(self.right_frame, text="Glycosylate")
       self.glycosylateAll = tk.Button(self.right_frame, text="Glycosylate all")
       self.clashes = tk.Button(self.right_frame, text="Remove clashes")
       # Layout right frame       
       i = 0
       self.w_logo.grid(column = 0, row =  i); i+=1
       self.chain_label.grid(column = 0, row =  i); i+=1
#       self.chain_menu.grid(column = 0, row =  i); i+=1
       self.sequon_label.grid(column = 0, row =  i); i+=1
#       self.sequon_menu.grid(column = 0, row =  i); i+=1
       self.glycan_label.grid(column = 0, row =  i); i+=1
       self.glycan_2D.grid(column = 0, row =  i); i+=1
       self.glycosylate.grid(column = 0, row =  i); i+=1
       self.glycosylateAll.grid(column = 0, row =  i); i+=1
       
            



       
       
if __name__ == "__main__":
    glycogui = GlycosylatorGUI()
    glycogui.mainloop()
       