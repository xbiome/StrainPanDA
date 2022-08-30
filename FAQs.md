# StrainPanDA - A strain analysis pipeline based on pangenome
# FAQs (Frequently asked questions)


## Errors pulling docker images


### Reason1: Did not login in docker.

Register an account in [docker hub](https://hub.docker.com/).

Login as the following [instruction](https://www.runoob.com/docker/docker-login-command.html)

<br>

### Reason2: the authorization of docker is not set correctly to the user.

Solution1: sudo docker pull

```
sudo docker pull yuxiangtan/strainpanda-mapping:dev
sudo docker tag yuxiangtan/strainpanda-mapping:dev strainpanda-mapping:dev
```

Solution2: add docker group for users by root. See [here](https://docs.docker.com/engine/install/linux-postinstall/) for details.

### 原因三：国内源太差

可以如图设置源路径，能较大程度提升下载速度。

## Failing to run PanPhlAN

### Main reason: the name of reference database is not correct.

The name of referece database must following the correct rule.

To double check, you can `ls -l <the path of error log>`

For example as the following:



## Setup the number of processors.


You can set it by changing the cpus value in the conf/base.config file, following [this](https://www.nextflow.io/docs/latest/process.html#cpus )


线程参数暂时没有开放给用户使用，默认使用8个线程
感兴趣的可以根据nextflow文档自行修改conf/base.config 中cpus的预设值 
https://www.nextflow.io/docs/latest/process.html#cpus 


## Get the functional annotation of gene families 

You can get the annotation of gene families in the prebuild databases by using the "anno.csv" files. In the real practice, you can use "merge" function in R to get the annotation of the gene family profile.
