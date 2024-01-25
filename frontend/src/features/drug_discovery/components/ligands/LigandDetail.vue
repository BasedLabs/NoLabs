<template>
     <q-card flat bordered v-if="!ligand.data">
        <q-card-section>
            <div class="text-caption">Smiles String</div>
        </q-card-section>
        <q-separator />

        <q-card-section>
            <q-skeleton :type="type" />
        </q-card-section>
        <q-card-section>
            <div class="text-caption">3D View</div>
        </q-card-section>
        <q-skeleton height="200px" square />
        <q-separator />
    </q-card>
    <q-card v-if="ligand.data">
        <q-card-section>
            <div class="text-h6">Ligand: {{ ligand.metaData.ligandName }}</div>
            <div ref="viewerContainer" class="ligand-viewer"></div>
        </q-card-section>
    </q-card>
</template>
  
<script>
import { QCard, QCardSection } from 'quasar';
import * as NGL from 'ngl';

export default {
    name: 'LigandDetail',
    components: {
        QCard,
        QCardSection,
    },
    props: {
        ligandMetadata: {
            ligand_id: String,
            ligandName: String,
        },
    },
    data() {
        return {
            viewer: null,
            ligand: {
                metaData: this.ligandMetadata,
                data: null
            }
        };
    },
    methods: {
        createViewer() {
            if (this.viewer || !this.ligand.ligandSdf) {
                return;
            }

            this.viewer = new NGL.Stage(this.$refs.viewerContainer);
            const stringBlob = new Blob([this.ligand.ligandSdf], { type: 'text/plain' });

            this.viewer.loadFile(stringBlob, { ext: 'sdf' }).then((component) => {
                component.addRepresentation('ball+stick');
                component.autoView();
            });
        },
    },
    mounted() {
        this.createViewer();
    },
    beforeUnmount() {
        if (this.viewer) {
            this.viewer.dispose();
        }
    },
};
</script>
  
<style>
.ligand-detail-page {
    /* Layout and spacing adjustments */
    margin-top: 20px;
}

.ligand-viewer {
    width: 100%;
    height: 400px;
    /* Adjust as needed */
}
</style>
  