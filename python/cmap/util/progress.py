'''
classes for determinate and indeterminate update bars
Created on Apr 10, 2012

@author: cflynn, wangz10
'''
import sys
from platform import platform
IS_WINDOWS = True
if not platform(terse=1).startswith('Windows'):
    import blessings
    IS_WINDOWS = False

from threading import Timer

class DeterminateProgressBar(object):
    '''
    provides an interface for determinate progress bars
    '''


    def __init__(self,name):
        '''
        Constructor
        '''
        self.name = name
        if not IS_WINDOWS:
            self.term = blessings.Terminal()
    
    def update(self,message, progress, total):
        '''
        update the update displayed on screen
        '''
        percent = float(progress)/total*100
        if not IS_WINDOWS:
            name_string = self.term.yellow(self.name)
            sys.stdout.write('\r' + name_string + ':%s  [%s] %.2f%%' %(
                                                     message,
                                                      '#'*(int(round(percent/10))), 
                                                      percent))
        else:
            sys.stdout.write('\r' + ':%s  [%s] %.2f%%' %(
                                                     message,
                                                      '#'*(int(round(percent/10))), 
                                                      percent))

        sys.stdout.flush()
        
    def show_message(self,message):
        '''
        displays the current message on screen until cleared by another class method
        '''
        self.clear()
        if not IS_WINDOWS:
            name_string = self.term.yellow(self.name)
            sys.stdout.write('\r' + name_string + ':%s' %(message,))
        else:
            sys.stdout.write('\r' + ':%s' %(message,))
        sys.stdout.flush()
    
    def clear(self):
        '''
        clears the screen
        '''
        try:
            if not IS_WINDOWS:
                sys.stdout.write('\r' +  ' ' * self.term.width)
            else:
                sys.stdout.write('\r' +  ' ' )
        except TypeError:
            sys.stdout.write('\r' +  ' ' * 1000)
        sys.stdout.flush()

class IndeteriminateProgressBar(object):
    '''
    provides an interface for indeterminate progress bars
    '''
    
    def __init__(self,name):
        '''
        Constructor
        '''
        self.name = name
        self.on = False
    
    def start(self):
        '''
        start the indeterminate progress bar
        '''
        self.on = True
        self.animate(0)
        
    
    def animate(self,i):
        '''
        underlying animate function for the progress bar
        '''
        while self.on:
            sys.stdout.write( '\r' + ( '.' * i ) + '   ' )
            sys.stdout.flush()
            Timer( 2, self.animate, ( 0 if i == 3 else i + 1, ) ).start()
    
    def stop(self):
        '''
        stop the indeterminate progress bar
        '''
        self.on = False
    
    def test(self):
        '''
        test the progress bar
        '''
        self.start()
        self.animate(0)
        for i in range(10000): #@UnusedVariable
            pass
        self.stop()