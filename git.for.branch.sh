if [ "$3" == "" ];
then
	echo "
sh $0 
<git commit -m > 
<git add ,example: sh $0 'update readme' 'README.md */*'>
<branch name: master or gh-pages>
"
	exit
fi
set -vex
commit=$1
add=$2
branch=$3

git checkout $branch
git add $add
git commit -m "$commit"
git push -u origin $branch

exit

## for first use git
#git init
#git add *
#git commit -m 'update'
#git remote add origin https://github.com/orangeSi/ClustersPloter.git
#git push -u origin master
#git pull origin master

# config password for git
#git config credential.helper store, store mean forver remember password. git config credential.helper 'cache –timeout=3600' for temperate remember password
# git push orgin master 
# ......
# after this will remember the password

