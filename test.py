import Tkinter as tk

class Root(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.label = tk.Label(self, text="Hello World", padx=5, pady=5) 
        self.label.pack()
        
if __name__ == "__main__":
    root = Root()
    root.mainloop()