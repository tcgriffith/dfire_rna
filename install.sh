
echo "## add environment variable DFIRE_RNA_HOME into ~/.bashrc"
echo "## please restart the terminal"
DFIRE_RNA_HOME=`pwd`
user_home=`echo ~`
envsetnot=`env | awk -F= '{print $1}' | grep "DFIRE_RNA_HOME"`
if [ -z "$envsetnot" ];then
        cat   << EOF >> ${user_home}/.bashrc
export DFIRE_RNA_HOME=${DFIRE_RNA_HOME}
EOF

fi


