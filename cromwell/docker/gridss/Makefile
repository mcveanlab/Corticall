VERSION = 0.1.0
TAG = quay.io/corticall/gridss:$(VERSION)

all: build push

build:
	docker build -t $(TAG) .

push:
	docker push $(TAG)
