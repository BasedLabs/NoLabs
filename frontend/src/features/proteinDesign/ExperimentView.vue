<template>
    <div>
        <q-separator></q-separator>
        <q-layout container style="height: 100vh" class="shadow-2 rounded-borders">
            <q-header :class="$q.dark.isActive ? 'bg-secondary' : 'bg-black'">
                <q-toolbar>
                    <q-toolbar-title>Binder design experiment: {{ experimentName }} <q-btn round
                            @click="changeExperimentName" color="info" size="sm" icon="edit" /></q-toolbar-title>
                    <q-btn class="q-mx-md" color="positive" @click="downloadPdbFile()"
                        :disable="store.experiment == null || store.experiment?.pdbsContent.length === 0"
                        label="Download pdb file" />
                    <q-btn color="positive" size="md" label="Binder parameters" @click="parameters = !parameters" />
                </q-toolbar>
            </q-header>
            <q-page-container>
                <div class="row" v-if="store.experiment != null && store.experiment.pdbsContent.length > 0">
                    <div class="col q-ma-sm">
                        <q-card flat bordered class="my-card">
                            <q-card-section>
                                <div class="text-h6">{{ originProteinName }}</div>
                            </q-card-section>

                            <q-separator inset />

                            <q-card-section>
                                <div id="viewport1" style="width: 100%; min-height: 400px;"></div>
                            </q-card-section>
                        </q-card>
                    </div>
                    <div class="col q-ma-sm">
                        <q-card flat bordered class="my-card">
                            <q-card-section>
                                <div class="text-h6">Binder</div>
                            </q-card-section>

                            <q-separator inset />

                            <q-card-section>
                                <div id="viewport2" style="width: 100%; min-height: 400px;"></div>
                            </q-card-section>
                        </q-card>
                    </div>
                </div>

            </q-page-container>
        </q-layout>
        <q-dialog v-model="parameters" position="right" :maximized="true">
            <q-card>
                <q-card-section class="row items-center q-pb-none">
                    <div class="text-h6">Protein binder design parameters</div>
                    <q-space />
                    <q-btn icon="close" flat round dense v-close-popup />
                </q-card-section>
                <q-card-section>
                    <q-form @submit="onSubmit" class="q-gutter-md">
                        <q-input filled v-model="formData.contig" label="Contig" lazy-rules
                            :rules="[val => val && val.length > 0 || 'Please type something']">
                            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
                                Confirm your chains have the residue numbers you're looking to diffuse over. 100 - Diffuses
                                a monomer 100 residues long. 50-100 - Diffuses a hetero-oligomer of lengths 50 and 100.
                                5-15/A10-25/30-40 - Builds 5-15 residues N-terminally of A10-25 from the input pdb, followed
                                by 30-40 residues to its C-terminus. B1-100/0 100-100 - Generates 100 residue long binders
                                to residues 1-100 of chain B.
                            </q-tooltip>
                        </q-input>
                        <q-input filled v-model="formData.hotspots" label="Hotspots">
                            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
                                The model optionally readily learns that it should be making an interface which involving
                                these hotspot residues. Input is ChainResidueNumber: A100 for residue 100 on chain A.
                            </q-tooltip>
                        </q-input>
                        <q-input filled type="number" v-model="formData.numberOfDesigns" label="Number of desings"
                            lazy-rules :rules="[val => val && val > 0 || 'Please type something']">
                            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
                                Number of designs to generate
                            </q-tooltip>
                        </q-input>
                        <q-input filled type="number" v-model="formData.timesteps" label="Timesteps" lazy-rules
                            :rules="[val => val && val > 0 || 'Please type something']">
                            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
                                Desired iterations to generate structure.
                            </q-tooltip>
                        </q-input>
                        <q-file filled bottom-slots accept=".pdb" :rules="[val => val]" v-model="formData.pdbFile"
                            label=".pdb file" counter>
                            <template v-slot:prepend>
                                <q-icon name="cloud_upload" @click.stop.prevent />
                            </template>
                            <template v-slot:append>
                                <q-icon name="close" @click.stop.prevent="formData.pdbFile = null" class="cursor-pointer" />
                            </template>

                            <template v-slot:hint>
                                Upload .pdb file with protein structure
                            </template>
                        </q-file>
                        <div>
                            <q-btn label="Submit" size="large" type="submit" color="positive" />
                        </div>
                    </q-form>
                </q-card-section>
            </q-card>
        </q-dialog>
    </div>
</template>
  
<script lang="ts">
import { defineComponent, ref } from 'vue'
import useProteinDesignStore from './storage';
import { QVueGlobals, useQuasar, QSpinnerOrbit } from 'quasar';

export default defineComponent({
    name: 'ProteinDesignExperimentView',
    data() {
        const store = useProteinDesignStore();
        const views = {
            default: { key: 'default', title: 'Default representation' },
            cartoon: { key: 'cartoon', title: 'Cartoon' },
            backbone: { key: 'backbone', title: 'Backbone' },
            ballsAndSticks: { key: 'ball+stick', title: 'Balls and sticks' },
            contact: { key: 'contact', title: 'Contact' },
            helixorient: { key: 'helixorient', title: 'Helixorient' },
            hyperball: { key: 'hyperball', title: 'Hyperball' },
            licorice: { key: 'licorice', title: 'Licorice' },
            ribbon: { key: 'ribbon', title: 'Ribbon' },
            rope: { key: 'rope', title: 'Rope' },
            surface: { key: 'surface', title: 'Surface' },
            spacefill: { key: 'spacefill', title: 'Spacefill' },
            unitcell: { key: 'unitcell', title: 'Unitcell' }
        };

        return {
            store,
            parameters: false,
            formData: {
                pdbFile: null,
                contig: '',
                numberOfDesigns: 1,
                timesteps: 50,
                hotspots: ''
            },
            experimentName: store.experiment === null ? 'Experiment name' : store.experiment!.name,
            quasar: null as unknown as QVueGlobals,
            views: views,
            viewsItems: Object.keys(views),
            selectedView: null,
            selectedTableIndex: 0,
            selectedPdb: null,
            stage: null,
            pdbComponent: null,
            ligandComponent: null,
            originProteinName: '',
        }
    },
    methods: {
        async onSubmit() {
            this.quasar.loading.show({
                spinner: QSpinnerOrbit,
                message: 'Running AI models. This can take a couple of minutes'
            });

            this.parameters = !this.parameters;
            await this.store.inference({
                experimentId: this.store.experiment?.id,
                experimentName: this.store.experiment?.name as string,
                pdbFile: this.formData.pdbFile as unknown as Blob,
                contig: this.formData.contig,
                numberOfDesigns: this.formData.numberOfDesigns,
                timesteps: this.formData.timesteps,
                hotspots: this.formData.hotspots
            })
            var read = new FileReader();

            read.readAsBinaryString(this.formData.pdbFile! as File);

            read.onloadend = () => {
                this.render('viewport1', read.result as string);
            }

            if (this.store.experiment != null) {
                this.render('viewport2', this.store.experiment!.pdbsContent[0]);
            }

            this.originProteinName = (this.formData.pdbFile! as File).name;
            this.quasar.loading.hide();
        },
        setView(viewKey: string) {
            this.selectedView = this.views[viewKey];
            if (this.selectedView === this.views.default) {
                this.render('viewport2',);
                return;
            }
            this.pdbComponent.removeAllRepresentations();
            this.pdbComponent.addRepresentation(this.selectedView.key);
        },
        render(viewer: string, pdbContent: string) {
            setTimeout(() => {
                // Render 3d structure
                document.getElementById(viewer).innerHTML = '';
                const stage = new NGL.Stage(viewer);
                stage.setParameters({ backgroundColor: 'black' });
                const proteinFileContentBlob = new Blob([pdbContent], { type: 'text/plain' });
                const proteinFile = new File([proteinFileContentBlob], 'protein.pdb', { type: 'text/plain' });
                stage.loadFile(proteinFile, { defaultRepresentation: true }).then((component) => {
                });
            }, 200);
        },
        downloadPdbFile() {
            const filename = `${this.store.experiment.name}.pdb`;
            const blob = new Blob([this.store.experiment?.pdbsContent[0]], { type: 'text/plain' });
            if (window.navigator.msSaveOrOpenBlob) {
                window.navigator.msSaveBlob(blob, filename);
            } else {
                const elem = window.document.createElement('a');
                elem.href = window.URL.createObjectURL(blob);
                elem.download = filename;
                document.body.appendChild(elem);
                elem.click();
                document.body.removeChild(elem);
            }
        },
        changeExperimentName() {
            this.quasar.dialog({
                title: 'Prompt',
                message: 'Enter new experiment name',
                prompt: {
                    model: '',
                    type: 'text' // optional
                },
                cancel: true,
                persistent: true
            }).onOk(async data => {
                await this.store.changeExperimentName(this.store.experiment?.id as string, data);
                this.experimentName = data;
            }).onCancel(() => {
                // console.log('>>>> Cancel')
            }).onDismiss(() => {
                // console.log('I am triggered on both OK and Cancel')
            })
        }
    },
    async mounted() {
        const experimentId = this.$route.params.experimentId as string;
        this.quasar = useQuasar();
        this.quasar.loading.show({
            spinner: QSpinnerOrbit,
            message: `Experiment ${experimentId}`
        })
        await this.store.getExperiment(experimentId);
        setTimeout(() => {
            this.quasar.loading.hide()
        }, 500);

        if (this.store.experiment && this.store.experiment.pdbsContent.length === 0) {
            this.parameters = true;
        }
    }
})
</script>
  