from Tkinter import * 

def onObjectClick(event):                  
        print 'Clicked', event.x, event.y, event.widget
        item = event.widget.find_closest(event.x, event.y)   
        print event.widget.gettags(item)

root = Tk()
canv = Canvas(root, width=100, height=100)
obj1 = canv.create_text(50, 30, text='Click me one', tags = 0)
obj2 = canv.create_text(50, 70, text='Click me two', tags = 1)

canv.tag_bind(obj1, '<Double-1>', onObjectClick)        
canv.tag_bind(obj2, '<Double-1>', onObjectClick)        
canv.pack()
root.mainloop()
