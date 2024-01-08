<script>
export default {
    props: ['experiment'],
    data() {
        const formatText = (prob) => {
            return prob.toFixed(2) * 100;
        }
        const localisation = this.experiment.data.localisation;
        const highlightData = {
            mithochondria: {
                image: 'localisationImages/mithochondria.png',
                controlElementId: 'mithochondria-list-item',
                text: `Mithochondria ${formatText(localisation.mithochondria)}%`
            },
            nucleus: {
                image: 'localisationImages/nucleus.png',
                controlElementId: 'nucleus-list-item',
                text: `Nucleus ${formatText(localisation.nucleus)}%`
            },
            cytoplasm: {
                image: 'localisationImages/cytoplasm.png',
                controlElementId: 'cytoplasm-list-item',
                text: `Cytoplasm ${formatText(localisation.cytoplasm)}%`
            },
            other: {
                image: 'localisationImages/original.png',
                controlElementId: 'other-proteins-item',
                text: `Other proteins ${formatText(localisation.other)}%`
            },
            extracellular: {
                image: 'localisationImages/original.png',
                controlElementId: 'extracellular-proteins-item',
                text: `Extracellular ${formatText(localisation.extracellular)}%`
            },
            original: {
                image: 'localisationImages/original.png'
            }
        };
        const obj = {
            highlightData: highlightData,
            activeObject: highlightData.original
        }
        return obj;
    },
    methods: {
        onMouseOver(el) {
            this.activeObject = this.highlightData.original;

            let anyActive = false;
            for (const key in this.highlightData) {
                const data = this.highlightData[key];
                // If we are over the list item, render corresponding image
                if (data.controlElementId == el.target.id) {
                    this.activeObject = data;
                    anyActive = true;
                }
            }

            // If not any active - render default
            if (!anyActive) {
                this.activeObject = this.highlightData.original;
            }
        },
        onMouseLeave(el) {
            this.activeObject = this.highlightData.original;
        }
    }
}
</script>

<template>
    <div class="col-md-12" id="localisation">
        <img id="img" :src="activeObject.image" class="localisation-image">
        <p>Hover on list</p>
        <ul class="list-group localisation-probs ">
            <li id="mithochondria-list-item" class="list-group-item"
                :class="highlightData.mithochondria === activeObject ? 'active' : ''" v-on:mouseover="onMouseOver"
                v-on:mouseleave="onMouseLeave">{{ highlightData.mithochondria.text }}</li>
            <li id="nucleus-list-item" class="list-group-item"
                :class="highlightData.nucleus === activeObject ? 'active' : ''" v-on:mouseover="onMouseOver"
                v-on:mouseleave="onMouseLeave">{{ highlightData.nucleus.text }}</li>
            <li id="cytoplasm-list-item" class="list-group-item"
                :class="highlightData.cytoplasm === activeObject ? 'active' : ''" v-on:mouseover="onMouseOver"
                v-on:mouseleave="onMouseLeave">{{ highlightData.cytoplasm.text }}</li>
            <li id="other-proteins-item" class="list-group-item"
                :class="highlightData.other === activeObject ? 'active' : ''" v-on:mouseover="onMouseOver"
                v-on:mouseleave="onMouseLeave">{{ highlightData.other.text }}</li>
            <li id="extracellular-proteins-item" class="list-group-item"
                :class="highlightData.extracellular === activeObject ? 'active' : ''" v-on:mouseover="onMouseOver"
                v-on:mouseleave="onMouseLeave">{{ highlightData.extracellular.text }}</li>
        </ul>
    </div>
</template>
