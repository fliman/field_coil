CC = /opt/local/bin/x86_64-apple-darwin14-g++-mp-4.9



run:main.o
	${CC} $< -o $@
