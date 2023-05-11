#!/bin/bash

ALL="../res/BOOST/ME/"

SIZE400="${ALL}400/"
SIZE800="${ALL}800/"
SIZE1200="${ALL}1200/"
SIZE1600="${ALL}1600/"

MODEL1="${SIZE400}1/ ${SIZE400}2/ ${SIZE400}3/ ${SIZE800}1/ ${SIZE800}2/ ${SIZE800}3/ ${SIZE1200}1/ ${SIZE1200}2/ ${SIZE1200}3/ ${SIZE1600}1/ ${SIZE1600}2/ ${SIZE1600}3/"
MODEL2="${SIZE400}4/ ${SIZE400}5/ ${SIZE400}6/ ${SIZE800}4/ ${SIZE800}5/ ${SIZE800}6/ ${SIZE1200}4/ ${SIZE1200}5/ ${SIZE1200}6/ ${SIZE1600}4/ ${SIZE1600}5/ ${SIZE1600}6/"
MODEL3="${SIZE400}7/ ${SIZE400}8/ ${SIZE400}9/ ${SIZE800}7/ ${SIZE800}8/ ${SIZE800}9/ ${SIZE1200}7/ ${SIZE1200}8/ ${SIZE1200}9/ ${SIZE1600}7/ ${SIZE1600}8/ ${SIZE1600}9/"
MODEL4="${SIZE400}10/ ${SIZE400}11/ ${SIZE400}12/ ${SIZE800}10/ ${SIZE800}11/ ${SIZE800}12/ ${SIZE1200}10/ ${SIZE1200}11/ ${SIZE1200}12/ ${SIZE1600}10/ ${SIZE1600}11/ ${SIZE1600}12/"

MAF01="${SIZE400}1/ ${SIZE400}4/ ${SIZE400}7/ ${SIZE400}10/ ${SIZE800}1/ ${SIZE800}4/ ${SIZE800}7/ ${SIZE800}10/ ${SIZE1200}1/ ${SIZE1200}4/ ${SIZE1200}7/ ${SIZE1200}10/ ${SIZE1600}1/ ${SIZE1600}4/ ${SIZE1600}7/ ${SIZE1600}10/"
MAF02="${SIZE400}2/ ${SIZE400}5/ ${SIZE400}8/ ${SIZE400}11/ ${SIZE800}2/ ${SIZE800}5/ ${SIZE800}8/ ${SIZE800}11/ ${SIZE1200}2/ ${SIZE1200}5/ ${SIZE1200}8/ ${SIZE1200}11/ ${SIZE1600}2/ ${SIZE1600}5/ ${SIZE1600}8/ ${SIZE1600}11/"
MAF04="${SIZE400}3/ ${SIZE400}6/ ${SIZE400}9/ ${SIZE400}12/ ${SIZE800}3/ ${SIZE800}6/ ${SIZE800}9/ ${SIZE800}12/ ${SIZE1200}3/ ${SIZE1200}6/ ${SIZE1200}9/ ${SIZE1200}12/ ${SIZE1600}3/ ${SIZE1600}6/ ${SIZE1600}9/ ${SIZE1600}12/"

python aggregate_results.py --csv-directories ${ALL} --save-as ../res/BOOST/ALL.csv --recursive

python aggregate_results.py --csv-directories ${SIZE400} --save-as ../res/BOOST/SIZE_400.csv --recursive
python aggregate_results.py --csv-directories ${SIZE800} --save-as ../res/BOOST/SIZE_800.csv --recursive
python aggregate_results.py --csv-directories ${SIZE1200} --save-as ../res/BOOST/SIZE_1200.csv --recursive
python aggregate_results.py --csv-directories ${SIZE1600} --save-as ../res/BOOST/SIZE_1600.csv --recursive

python aggregate_results.py --csv-directories ${MODEL1} --save-as ../res/BOOST/MODEL_1.csv
python aggregate_results.py --csv-directories ${MODEL2} --save-as ../res/BOOST/MODEL_2.csv
python aggregate_results.py --csv-directories ${MODEL3} --save-as ../res/BOOST/MODEL_3.csv
python aggregate_results.py --csv-directories ${MODEL4} --save-as ../res/BOOST/MODEL_4.csv

python aggregate_results.py --csv-directories ${MAF01} --save-as ../res/BOOST/MAF_01.csv
python aggregate_results.py --csv-directories ${MAF02} --save-as ../res/BOOST/MAF_02.csv
python aggregate_results.py --csv-directories ${MAF04} --save-as ../res/BOOST/MAF_04.csv