APPS = \
	send-receive trapezoid mv-mult

send-receive: send-receive.cpp
	mpic++ $? -o $@

trapezoid: trapezoid.cpp
	mpic++ $? -o $@

mv-mult: mv-mult.cpp
	mpic++ $? -o $@


clean:
	rm $(APPS)

all: $(APPS)