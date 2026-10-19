if [ ! -f subjList.txt ]; then
	ls | grep ^sub- > subjList.txt
fi
