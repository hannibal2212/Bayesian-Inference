library(igraph)


query<-c() # n elements
startnode<-1   # the node to start
query[1]<-startnode;
distance<-rep(Inf,nrow(givengraph))
distance[1]<-0;
writepointer<-2 # points the element in query to write 
readpointer<-1  # points the element in query to read 
visited <- c()

while(writepointer!= readpointer){
  node<-query[readpointer]
  for(j in which(!(seq_along(givengraph[node,]) %in% query))){
    if( givengraph[node,j]==1){
      query[writepointer]<-j;
      distance[writepointer]<-distance[readpointer]+1;
      writepointer = writepointer + 1
    }
    
    
  }
  
  readpointer = readpointer + 1
}


adjm1<-matrix(sample(0:1,100,replace=TRUE,prob=c(0.9,01)),nc=10)
givengraph<-adjm1

givengraph[3,1] == 0



givengraph<-matrix(c(0,1,0,0,0,0,
         1,0,1,0,0,0,
         0,1,0,1,1,0,
         0,0,1,0,1,0,
         0,0,1,1,0,1,
         0,0,0,0,1,0),nrow = 6, ncol = 6, byrow = T)



givengraph<-matrix(c(0,	1,	1,	1,	0,	0,	0,	0,	0,	0,
  1,	0,	0,	0,	1,	1,	1,	0,	0,	0,
  1,	0,	0,	0,	1,	0,	0,	0,	0,	0,
  1,	0,	0,	0,	0,	1,	0,	0,	0,	0,
  0,	1,	1,	0,	0,	0,	0,	0,	1,	1,
  0,	1,	0,	1,	0,	0,	0,	0,	1,	1,
  0,	1,	0,	0,	0,	0,	0,	1,	0,	0,
  0,	0,	0,	0,	0,	0,	1,	0,	0,	0,
  0,	0,	0,	0,	1,	1,	0,	0,	0,	0,
  0,	0,	0,	0,	1,	1,	0,	0,	0,	0), nrow = 10, ncol = 10 , byrow = T)


which(givengraph[5,]==1)


setdiff(givengraph[5,],givengraph[5,2])


givengraph[5,-c(2)]


dist1<-rep(Inf,n)
dist2<-rep(Inf,n)
dist2[r]<-dist1[r]<-0
v<-1:n
vt<-c(r)
while(length(v)!=length(vt)){
  for(i in length(vt)){
    for (j in which(!(v %in% vt))){
      dist2[j]<-gg[vt[i],j]
      
      
      vt<-c(vt,k)
    }
    dist1[k]<-min(dist2)
    k<-which(dist2 == min(dist2))
    if(dist2[k]<dist1[k]){dist1[k]<-dist2[k]}
  }
}


dist_first<-rep(Inf,6)
dist_subs<-rep(Inf,6)
dist_subs[2]<-0
vt<-c(2)
v<-1:6
while(length(vt)!=length(v)){
  for (i in 1:length(vt)){
    dist_subs[which(!(v %in% vt))]<-gg[vt[i],-vt]
    for(j in v){
      if(dist_subs[j]<dist_first[j]){
        dist_first[j]<-dist_subs[j]
        
      }
      
    }
    vt<-c(vt,which (dist_subs == min(dist_subs[-vt]))) 
    
  }
}


gg<-matrix(c(0,1,3,Inf,Inf,3,
             1,0,5,1,Inf,Inf,
             3,5,0,2,1,Inf,
             Inf,1,2,0,4,Inf,
             Inf,Inf,1,4,0,5,
             3,Inf,Inf,Inf,5,0),nrow = 6,ncol = 6,byrow = T)



if(dist_subs[j] == dist_first[j]){
  dist_first1<-dist_first
  dist_first[j]<-dist_subs[j]
  dist_first2<-dist_first
}








