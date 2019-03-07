
echo "## set DFIRE_RNA_HOME environment"
DFIRE_RNA_HOME=`pwd`
user_home=`echo ~`
envsetnot=`env | awk -F= '{print $1}' | grep "DFIRE_RNA_HOME"`
if [ -z "$envsetnot" ];then
        cat   << EOF >> ${user_home}/.bashrc
export DFIRE_RNA_HOME=${DFIRE_RNA_HOME}
EOF
        source ${user_home}/.bashrc
fi

make 
