VERSION = 0.1.1
TAG = quay.io/corticall/utils:$(VERSION)

all: build push

build:
	docker build -t $(TAG) .

push:
	docker push $(TAG)
