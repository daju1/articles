#!/usr/bin/env bash

set -e
set -x


if [ ! -d  ~/bin/ ]; then
	mkdir  ~/bin/
fi

cd git-trac-command
ln -s `pwd`/git-trac ~/bin/