var map = function(){
  var data = {'w':this.w, 'q2':this.q2 }
  emit(data, 1);
}
var reduce = function(key, values){
  var res = 0;
  values.forEach(function(v){ res += 1});
  return {count: res};
}

db.events.mapReduce(map, reduce, { out: 'wvq2', results: "wvq2" });
