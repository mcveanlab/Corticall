VERSION = 0.2.0
TAG = quay.io/corticall/mccortex:$(VERSION)

all: build push

build:
	docker build -t $(TAG) .

push:
	docker push $(TAG)
