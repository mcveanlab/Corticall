VERSION = 0.1.9
TAG = quay.io/corticall/corticall:$(VERSION)

all: build push

build:
	cp ../../../dist/corticall.jar .
	docker build -t $(TAG) .

push:
	docker push $(TAG)
