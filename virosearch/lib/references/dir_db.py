import os
import re
import pathlib


###############################################################################
#                                    DirDB                                    #
###############################################################################

# This class is responsible for setting up a structured ViroSearch database
# directory and automatically downloading a list of reference sequence files,
# given a list of accession IDs. If the ViroSearch directory already exists, we
# will be appending new reference files to the directory structure. As this
# class scans the directory structure, a new class instance should be invoked
# prior to making any decision---i.e. if files are added to the directory
# structure after being called, a preexisting object will not see them.


class DirDB:

    def __init__(self, dir_, verbose=True):

        self.dir_ = dir_
        self.verbose = verbose
        self.references = dict()
        self.dirs = set()

        # Recursively traverse the given ViroSearch directory and collect file
        # paths of interest. As the directory structure is standardized, we may
        # use this class to evaluate and navigate the file structure of the
        # directory. A generic example of this tree follows:

        # ref_EBV
        # └── referencesDB
        #     ├── NC_007605.1
        #     │   ├── NC_007605.1.cds.faa
        #     │   ├── NC_007605.1.cds.fna
        #     │   ├── NC_007605.1.fna
        #     │   ├── NC_007605.1.fna.amb
        #     │   ├── NC_007605.1.fna.ann
        #     │   ├── NC_007605.1.fna.bwt
        #     │   ├── NC_007605.1.fna.cmd
        #     │   ├── NC_007605.1.fna.fai
        #     │   ├── NC_007605.1.fna.pac
        #     │   ├── NC_007605.1.fna.sa
        #     │   └── NC_007605.1.gbk
        #     └── NC_009334.1
        #         ├── NC_009334.1.cds.faa
        #         ├── NC_009334.1.cds.fna
        #         ├── NC_009334.1.fna
        #         ├── NC_009334.1.fna.amb
        #         ├── NC_009334.1.fna.ann
        #         ├── NC_009334.1.fna.bwt
        #         ├── NC_009334.1.fna.cmd
        #         ├── NC_009334.1.fna.pac
        #         ├── NC_009334.1.fna.sa
        #         └── NC_009334.1.gbk

        # NOTE: The bwa index files (.amb, .ann, .bwt, .pac, .sa) may or may
        # not be present.

        dir_db = pathlib.Path(self.dir_).resolve()

        for dir_path, dir_names, file_names in os.walk(dir_db):

            if dir_names:
                for dir_ in dir_names:
                    self.dirs.add(dir_)

            if file_names:

                for file_ in file_names:

                    file_path = pathlib.Path(dir_path, file_)
                    file_name = file_path.name

                    accid = file_name
                    accid = re.sub('.cds.fna', '', accid)
                    accid = re.sub('.cds.faa', '', accid)
                    accid = re.sub('.gbk', '', accid)
                    accid = re.sub('.fna', '', accid)
                    accid = re.sub('.amb', '', accid)
                    accid = re.sub('.ann', '', accid)
                    accid = re.sub('.bwt', '', accid)
                    accid = re.sub('.pac', '', accid)
                    accid = re.sub('.sa', '', accid)

                    # Populate the instance's reference dictionary with file
                    # paths.

                    file_path = file_path.resolve().as_posix()

                    if accid not in self.references.keys():
                        if file_name.endswith('.gbk'):
                            self.references.update({accid: {'gbk': file_path}})
                        elif file_name.endswith('.cds.fna'):
                            self.references.update({accid: {'cds_fna$': file_path}})
                        elif file_name.endswith('.cds_faa'):
                            self.references.update({accid: {'cds_faa$': file_path}})
                        elif file_name.endswith('.fna'):
                            self.references.update({accid: {'fna': file_path}})
                        elif file_name.endswith('.amb'):
                            self.references.update({accid: {'bwa_amb': file_path}})
                        elif file_name.endswith('.ann'):
                            self.references.update({accid: {'bwa_ann': file_path}})
                        elif file_name.endswith('.bwt'):
                            self.references.update({accid: {'bwa_bwt': file_path}})
                        elif file_name.endswith('.pac'):
                            self.references.update({accid: {'bwa_pac': file_path}})
                        elif file_name.endswith('.sa'):
                            self.references.update({accid: {'bwa_sa': file_path}})
                    else:
                        if file_name.endswith('.gbk'):
                            self.references[accid]['gbk'] = file_path
                        elif file_name.endswith('.cds.fna'):
                            self.references[accid]['cds_fna'] = file_path
                        elif file_name.endswith('.cds.faa'):
                            self.references[accid]['cds_faa'] = file_path
                        elif file_name.endswith('.fna'):
                            self.references[accid]['fna'] = file_path
                        elif file_name.endswith('.amb'):
                            self.references[accid]['bwa_amb'] = file_path
                        elif file_name.endswith('.ann'):
                            self.references[accid]['bwa_ann'] = file_path
                        elif file_name.endswith('.bwt'):
                            self.references[accid]['bwa_bwt'] = file_path
                        elif file_name.endswith('.pac'):
                            self.references[accid]['bwa_pac'] = file_path
                        elif file_name.endswith('.sa'):
                            self.references[accid]['bwa_sa'] = file_path

        # Deal with unset bwa index files.

        for accid in self.references.keys():

            if 'bwa_amb' not in self.references[accid].keys():
                self.references[accid]['bwa_amb'] = None

            if 'bwa_ann' not in self.references[accid].keys():
                self.references[accid]['bwa_ann'] = None

            if 'bwa_bwt' not in self.references[accid].keys():
                self.references[accid]['bwa_bwt'] = None

            if 'bwa_pac' not in self.references[accid].keys():
                self.references[accid]['bwa_pac'] = None

            if 'bwa_sa' not in self.references[accid].keys():
                self.references[accid]['bwa_sa'] = None

        return

    def for_accids(self):
        return self.references.keys()

    def gbk_file_path(self, accid):
        return self.references[accid]['gbk']

    def cds_fna_file_path(self, accid):
        return self.references[accid]['cds_fna']

    def cds_faa_file_path(self, accid):
        return self.references[accid]['cds_faa']

    def fna_file_path(self, accid):
        return self.references[accid]['fna']

    def bwa_amb_file_path(self, accid):
        return self.references[accid]['bwa_amb']

    def bwa_ann_file_path(self, accid):
        return self.references[accid]['bwa_ann']

    def bwa_bwt_file_path(self, accid):
        return self.references[accid]['bwa_bwt']

    def bwa_pac_file_path(self, accid):
        return self.references[accid]['bwa_pac']

    def bwa_sa_file_path(self, accid):
        return self.references[accid]['bwa_sa']

    def dir_path(self):
        return self.dir_.resolve()

    def is_bwa_index(self, accid):

        # If any of the required bwa index files are missing, we will return
        # false.

        bool_ = True

        if self.bwa_amb_file_path(accid) is None:
            bool_ = False

        if self.bwa_ann_file_path(accid) is None:
            bool_ = False

        if self.bwa_bwt_file_path(accid) is None:
            bool_ = False

        if self.bwa_pac_file_path(accid) is None:
            bool_ = False

        if self.bwa_sa_file_path(accid) is None:
            bool_ = False

        return bool_

    def is_variants_db_dir(self):
        if 'variantsDB' in self.dirs:
            return True
        else:
            return False

    def is_references_db(self):
        if 'referencesDB' in self.dirs:
            return True
        else:
            return False

    def is_dir_path(self):
        return self.dir_.resolve().is_dir()


# __END__
