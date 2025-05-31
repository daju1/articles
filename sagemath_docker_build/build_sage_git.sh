#!/bin/sh

set -e
set -x

SAGE_VER=9.6.beta7
SAGE_DIR=sage

WORK_DIR=`pwd`

SAGE_INSTALL=/opt/sagemath-${SAGE_VER}

export SAGE_ROOT=$WORK_DIR/$SAGE_DIR
export SAGE_LOCAL=$SAGE_ROOT/local

TRAC_DIR=git-trac-command

if [ ! -d $TRAC_DIR ]; then
	git clone https://github.com/sagemath/git-trac-command.git
fi

#source git-trac-command/enable.sh


if [ ! -d $SAGE_DIR ]; then
	git clone https://github.com/sagemath/sage.git

	cd $SAGE_DIR

	# commit 78bfb6c7ad4e6c1e6df83d13c3789e2fb59a8d5d (tag: 9.6.beta7)

	git checkout 9.6.beta7


	# https://doc.sagemath.org/html/en/developer/manual_git.html
	git remote add trac https://github.com/sagemath/sagetrac-mirror.git -t master
	git remote set-url --push trac git@trac.sagemath.org:sage.git



	#git branch
	#  develop
	#  master
	#  t/32160/integrate_of_function_which_contains__x___floor_x__
	#  t/32161/derivative_of_the_symbolic_sum_of_function_of_two_variables
	#* t/32184/derivative_of_symbolic_integration_with_infinity_limit_fails

	git-trac checkout 32184
	git-trac checkout 32161
	git-trac checkout 32160

	# u/mkoeppe/boolean_symbolic_expressions
	git-trac checkout 31911

	git merge t/32184/derivative_of_symbolic_integration_with_infinity_limit_fails
	git merge t/32161/derivative_of_the_symbolic_sum_of_function_of_two_variables
	#git merge t/32160/integrate_of_function_which_contains__x___floor_x__
	git merge t/32160/public/ticket/32160

	cd ..
fi

cd $SAGE_DIR

./bootstrap

./configure  --prefix=$SAGE_LOCAL --disable-doc #--enable-build-as-root
# MAKE='make -jNUM' make
make -j$(nproc)
sudo make install

cd sage
./sage -i fricas


# https://stackoverflow.com/questions/40607370/make-pdf-in-jupyter-notebook-with-russian-letters
# $HOME/.local/share/jupyter/nbconvert/templates/latex/base.tex.j2

# \usepackage[T2A]{fontenc}
# \usepackage[utf8]{inputenc}
# \usepackage[english, russian]{babel}

# ((( super() )))
# \usepackage{wrapfig}
# \usepackage[rightcaption]{sidecap}
# \providecommand{\keywords}[1]{\textbf{\textit{Keywords:}} #1}

# \author{A. Yu. Drozdov}

