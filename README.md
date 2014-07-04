catalyzer
=========

```
           _        _                    
  ___ __ _| |_ __ _| |_   _ _______ _ __ 
 / __/ _` | __/ _` | | | | |_  / _ \ '__|
| (_| (_| | || (_| | | |_| |/ /  __/ |   
 \___\__,_|\__\__,_|_|\__, /___\___|_|   
     / ___________________/              
     \ \
      \ \
       \ \
        | |
        | |
        |_| 
```      

The catalyzer is a project to convert audio from the ultrasonic frequencies
down to the audible range in real time.  The theory is to capture data at
200khz (requires a special mic and pre-amp); FFT the data; divide all
frequencies by 4 and re-bin them; and finally inverse-FFT (and normalize) the
result. 

cmd:  arecord --device=hw:CARD=X20,DEV=0 -t raw -f S32_LE -c 2 -r 192000 | ./catalyze | aplay --device=hw:CARD=X20,DEV=0 -f S32_LE -c 2 -r 192000

To build:
	1) install fftw3 (in debian/ubuntu: apt-get install libfftw3-dev
	2) cd src; make
        3) be happy
