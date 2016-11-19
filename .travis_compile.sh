#echo "Download SeqAn..."
#curl http://packages.seqan.de/seqan-library/seqan-library-2.1.0.tar.xz | tar -xJf -
#echo "Copy SeqAn..."
#sudo cp -rf seqan-library-2.1.0/include/seqan /usr/local/include/
#echo "Done."
#
#echo "Patch include dir for Ubuntu..."
#sed -i 's|/usr/include|/usr/local/include|' setup.py
#echo "Done."

#python setup.py build
#python setup.py install --install-lib test
#nosetests test/simple_test.py

bash build_it
python test/test.py
