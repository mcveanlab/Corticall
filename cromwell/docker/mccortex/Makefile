VERSION = 0.3.0
TAG = quay.io/corticall/mccortex:$(VERSION)

all: build push

build:
	docker build -t $(TAG) .

push:
	docker push $(TAG)
