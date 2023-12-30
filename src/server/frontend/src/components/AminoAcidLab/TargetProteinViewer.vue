<script>
export default {
    props: ['target'],
    data() {
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
            views: views,
            viewsItems: Object.keys(views),
            selectedView: views.default,
            stage: null,
            pdbComponent: null,
            ligandComponent: null,
            selectedResidues: new Set(),
            sequence: [],
            selectionMode: false
        }
    },
    methods: {
        setView(evt) {
            const viewKey = evt.target.value;
            this.selectedView = this.views[viewKey];
            if (this.selectedView === this.views.default) {
                this.render();
                return;
            }
            this.pdbComponent.removeAllRepresentations();
            this.pdbComponent.addRepresentation(this.selectedView.key);
        },
        cleanStage() {
            if (this.stage)
                this.stage.removeAllComponents();
        },
        render() {
            setTimeout(() => {
                // Render 3d structure
                document.getElementById('viewport').innerHTML = '';
                this.cleanStage();
                this.stage = new NGL.Stage("viewport");
                this.stage.setParameters({ backgroundColor: 'white',  tooltip: true });
                const proteinFileContentBlob = new Blob([this.target.pdb], { type: 'text/plain' });
                const proteinFile = new File([proteinFileContentBlob], 'protein.pdb', { type: 'text/plain' });
                this.stage.loadFile(proteinFile, { defaultRepresentation: true }).then((component) => {
                    this.pdbComponent = component;
                    this.extractSequenceFromComponent(component);
                });
                this.stage.signals.clicked.add(this.handleStageClick);
            }, 200);
        },
        handleStageClick(pickingProxy) {
            if (this.selectionMode && pickingProxy && pickingProxy.atom) {
                this.toggleResidueSelection(pickingProxy.atom.residueIndex);
            }
        },
        highlightSelectedResidues() {
            if (!this.pdbComponent) return;

            // Remove only the specific representation for selected residues
            this.pdbComponent.removeRepresentation(this.selectedRepresentation);

            // Add new representation for the selected residues
            const residueIndices = Array.from(this.selectedResidues).join(" OR ");
            if (residueIndices) {
                this.selectedRepresentation = this.pdbComponent.addRepresentation("ball+stick", {
                    sele: `(${residueIndices})`,
                    color: "blue"
                });
            }
        },
        toggleSelectionMode() {
            this.selectionMode = !this.selectionMode;
        },
        toggleResidueSelection(index) {
            if (this.selectedResidues.has(index)) {
                this.selectedResidues.delete(index);
            } else {
                this.selectedResidues.add(index);
            }
            this.highlightSelectedResidues();
        },
        sendSelectedResidues() {
            const selectedResidueArray = Array.from(this.selectedResidues);
            console.log("Sending selected residues to backend:", selectedResidueArray);
            // Replace with your API call to send data to the backend
        },
        extractSequenceFromComponent(component) {
            // Extracting amino acid sequence from the PDB component
            let sequence = [];
            component.structure.eachResidue(function (residue) {
                if (residue.isProtein()) {
                    sequence.push(residue.getResname1());
                }
            });

            // Update sequence data property
            this.sequence = sequence;
        },
        downloadPdbFile(evt) {
            const filename = this.target.metadata.name + '.pdb';
            const blob = new Blob([this.target.pdb], { type: 'text/plain' });
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
        }
    },
    computed: {
        hasPdb() {
            return this.target.pdb != null;
        },
        fastaSequence() {
            // Assuming 'fasta' is a string of amino acids
            return this.target.fasta ? this.target.fasta : '';
        },
    },
    mounted() {
        this.render();
    }
}
</script>

<template>
    <div class="container-fluid">
        <div class="row">
            <div class="col-lg-8" v-if="hasPdb">
                <div id="viewport" style="width: 100%; height: 500px;"></div>
                <div class="amino-acid-sequence" style="width: 100%; color: black; background-color: #f8f8f8; padding: 10px;">
                    <span v-for="(acid, index) in sequence" :key="index"
                          :class="{ 'selected': selectedResidues.has(index) }"
                          @click="toggleResidueSelection(index)">
                        {{ acid }}
                    </span>
                </div>
            </div>
            <div v-else>
                <div class="amino-acid-sequence">{{ fastaSequence }}</div>
                <button class="btn btn-primary" @click="uploadPdb">Upload corresponding PDB</button>
                <button class="btn btn-secondary" @click="predictStructure">Predict 3D structure</button>
            </div>
            <div class="col-lg-4" v-if="hasPdb">
                <h4>Options</h4>
                <select class="form-select form-select-md mb-3" aria-label="Select a representation" @change="setView($event)">
                    <option v-for="viewKey in viewsItems" :key="viewKey" :value="viewKey">
                        {{ views[viewKey].title }}
                    </option>
                </select>
                <button type="button" @click="downloadPdbFile()"
                    class="btn btn-primary padding-top-button-group download-pdb-button">Download .pdb
                </button>
                <button type="button" @click="toggleSelectionMode" class="btn btn-secondary" style="margin-top: 10px;">
                    Select binding pocket
                </button>
                <button v-if="this.selectionMode" type="button" @click="sendSelectedResidues" class="btn btn-success" style="margin-top: 10px;">
                    Confirm selection
                </button>
            </div>
        </div>
        <div class="row">
            <div class="col-12" style="margin-top: 10px;">
                <button class="btn btn-danger" @click="$emit('close')">Close</button>
            </div>
        </div>
    </div>
</template>

<style>
.amino-acid-sequence {
    overflow-x: auto;
    white-space: nowrap;
    font-family: monospace;
    margin-bottom: 10px;
    border: 1px solid #ddd; /* Optional: adds a border around the sequence */
    border-radius: 5px; /* Optional: rounds the corners of the border */
}

.amino-acid-sequence span {
    display: inline-block;
    padding: 2px 5px;
    cursor: pointer;
    border-radius: 3px;
}

.amino-acid-sequence span.selected {
    background-color: blue;
    color: white;
}
</style>
