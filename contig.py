

class Contig():
    
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

        self.homo_dup_depth = []
        self.homo_non_dup_depth = []
    
        self.homo_dup_kmers = []
        self.dnd_ratio = []

        self.duplicated = []

    def calculate_dnd_ratio(self):

        for pos in range(len(self.homo_dup_depth)):
            # no homozygous kmers in this position
            if self.homo_dup_depth[pos] == 0 and self.homo_non_dup_depth[pos] == 0:
                self.dnd_ratio.append(0.5) # TODO: find a better way to handle no data
            else:
                # ie. percent of homozygous kmers that are duplicated
                self.dnd_ratio.append(self.homo_dup_depth[pos] / (self.homo_dup_depth[pos] + self.homo_non_dup_depth[pos]))    
    
    def plot_dnd_ratio(self):

        def moving_average(data, window_size):
            return np.convolve(data, np.ones(window_size) / window_size, mode='valid')
       
        moving_ave = moving_average(self.dnd_ratio, 1000)
        # print(self.dnd_ratio)
        # print(moving_ave)
        pos = [i for i in range(0, len(moving_ave))]

        if not os.path.exists("results"):
            os.makedirs("results")
            
        fig = px.scatter(x=pos, y=moving_ave, labels={'x': 'Position', 'y': '% duplicated kmers'})
        fig.write_image(f'results/{self.name}_dnd_ratio.png')
        fig.write_html(f'results/{self.name}_dnd_ratio.html')


        # pos = [i for i in range(0, len(self.dnd_ratio))]

        # fig = px.scatter(x=pos, y=self.dnd_ratio, labels={'x': 'Position', 'y': '% duplicated kmers'})
        # fig.write_image(f'results/{self.name}_dnd_ratio.png')
        # fig.write_html(f'results/{self.name}_dnd_ratio.html')

    def get_kmers(self, bam):
        """
        get kmers from a bam file
        """
        logging.info(f"reading bam: {bam} for kmers to {self.name}")
        cmd = f"samtools view {bam} -@ 8 '{self.name}'"  
        logging.info(cmd)
        
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

        while True:
            line = proc.stdout.readline()
            if not line:
                break

            # print("test:", line.decode('UTF-8'))
            line = line.decode('UTF-8').strip().split()
            self.homo_dup_kmers.append(line[0])

    def get_duplicated_sequence(self):
        pass
        # if not self.duplicated:
        #     return None
        # else:
        #     # If completely duplicated
        #     if self.duplicated[1] - self.duplicated[0] == len(self.sequence):
        #         return f">{self.name}\n{self.sequence}\n"
            
        #     # If 5' duplicated
        #     elif 0 in self.duplicated[0]:
        #         return f">{self.name}\n{self.sequence[0:self.duplicated[0]]}\n"
            
        #     # If 3' duplicated
        #     else :
        #         return f">{self.name}\n{self.sequence[self.duplicated[1]:]}\n"

    def get_non_duplicated_sequence(self):

        print(self.duplicated)
        if not self.duplicated:
            return f">{self.name}\n{self.sequence}"
        else:
            # If completely duplicated
            for interval in self.duplicated:
                if interval[1] - interval[0] == len(self.sequence):
                    return ""
            
            # If 5' duplicated
            for interval in self.duplicated:
                if 0 in interval:
                    return f"{self.name}\n{self.sequence[interval[1]:]}\n"
           
            # If 3' duplicated
            for interval in self.duplicated:
                if inverval[1] == len(self.sequence):
                    return f"{self.name}\n{self.sequence[0:interval[0]]}\n"

    def __str__(self):
        return f"contig: {self.name}"
    
    def __repr__(self):
        return f"contig: {self.name}"