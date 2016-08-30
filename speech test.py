import pyttsx
engine = pyttsx.init()
# engine.say('Sally sells seashells by the seashore.')
# engine.say('The quick brown fox jumped over the lazy dog.')
# engine.say('How now brown cow')
# engine.say('Peter piper picked a peck of pickled peppers')
import winsound
Freq = 2500 # Set Frequency To 2500 Hertz
Dur = 1000 # Set Duration To 1000 ms == 1 second
winsound.Beep(Freq,Dur)
while(1):
    Freq = 2500 # Set Frequency To 2500 Hertz
    Dur = 500 # Set Duration To 1000 ms == 1 second
    winsound.Beep(Freq,Dur)
    Freq = 1500 # Set Frequency To 2500 Hertz
    Dur = 500 # Set Duration To 1000 ms == 1 second
    winsound.Beep(Freq,Dur)
    engine.say('emergency')
    #engine.say ('Oh dear')
    engine.say('self destruct initiated')
    engine.runAndWait()
