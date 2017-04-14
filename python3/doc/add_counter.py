import glob

counter_code = '''
<!-- Start of StatCounter Code -->
<script type="text/javascript">
var sc_project=5680199; 
var sc_invisible=1; 
var sc_partition=63; 
var sc_click_stat=1; 
var sc_security="6e579550"; 
</script>

<script type="text/javascript"
src="http://www.statcounter.com/counter/counter.js"></script><noscript><div
class="statcounter"><a title="counter for iweb"
href="http://www.statcounter.com/iweb/" target="_blank"><img
class="statcounter"
src="http://c.statcounter.com/5680199/0/6e579550/1/"
alt="counter for iweb" ></a></div></noscript>
<!-- End of StatCounter Code -->
'''

def process_file(fn):
   lines = open(fn).readlines()
   for i in range(len(lines)):
      if lines[i].count('statcounter') > 0:
         return
   for i in range(len(lines)):
      if lines[i].count('<\body>') > 0:
         break
   lines.insert(i, counter_code)
   f = open(fn, 'w')
   for l in lines:
      f.write(l)
   f.close()
      
for fn in glob.glob("_build/html/*html"):
   process_file(fn)      

   
