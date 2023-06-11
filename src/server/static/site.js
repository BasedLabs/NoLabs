var localisation = function (probabilities) {
    window.onmousemove = function (e) {
        //console.log("mouse location:", e.clientX, e.clientY)
    };
    if (!document.getElementById('img'))
        throw Error('Declare img tag');
    const _createImage = (src) => {
        const image = new Image(500, 500);
        image.src = src;
        return image;
    };
    const img = document.getElementById('img');
    const obj = {
        hoverControlHighlight: {
            onmouseover: (el) => {
                const jqueryEl = $(el.target);
                jqueryEl.addClass('active');
                const localisationObjectKey = jqueryEl.data('localisationObjectKey');
                const highlightElement = obj.highlightData[localisationObjectKey];
                img.src = highlightElement.image.src;
            },
            onmouseleave: (el) => {
                $(el.target).removeClass('active');
                img.src = obj.highlightData.original.image.src;
            },
        },
        highlightData: {
            mithochondria: {
                image: _createImage('./static/images/mithochondria.png'),
                controlElementId: 'mithochondria-list-item',
                text: 'Mithochondria'
            },
            nucleus: {
                image: _createImage('./static/images/nucleus.png'),
                controlElementId: 'nucleus-list-item',
                text: 'Nucleus'
            },
            cytoplasm: {
                image: _createImage('./static/images/cytoplasm.png'),
                controlElementId: 'cytoplasm-list-item',
                text: 'Cytoplasm'
            },
            other: {
                image: _createImage('./static/images/original.png'),
                controlElementId: 'other-proteins-item',
                text: 'Other proteins'
            },
            extracellular: {
                image: _createImage('./static/images/original.png'),
                controlElementId: 'extracellular-proteins-item',
                text: 'Extracellular proteins'
            },
            original: {
                image: _createImage('./static/images/original.png')
            },
        },
        render: function () {
            const highlightElements = [];
            for (const key in obj.highlightData) {
                const data = obj.highlightData[key];
                const elId = data.controlElementId;
                const el = $('#' + elId);
                if (probabilities[key] === undefined)
                    probabilities[key] = 0.000
                el.append(data.text + ' ' + (probabilities[key].toFixed(2) * 100) + '%');
                el.hover(obj.hoverControlHighlight.onmouseover,
                    obj.hoverControlHighlight.onmouseleave);
                el.data('localisationObjectKey', key);
                highlightElements.push(el);
            }

            img.src = obj.highlightData.original.image.src;
        }
    }
    return obj;
}
var folding = function (proteinFileContent) {
    const stage = new NGL.Stage("viewport");
    stage.setParameters({backgroundColor: 'white'})
    const obj = {
        views: {
            default: 'default',
            cartoon: 'cartoon',
            backbone: 'backbone',
            ballsAndSticks: 'ball+stick',
            contact: 'contact',
            helixorient: 'helixorient',
            hyperball: 'hyperball',
            licorice: 'licorice',
            ribbon: 'ribbon',
            rope: 'rope',
            surface: 'surface',
            spacefill: 'spacefill',
            unitcell: 'unitcell'
        },
        reload: () => {
            stage.removeAllComponents();
            const proteinFileContentBlob = new Blob([proteinFileContent], {type: 'text/plain'});
            const proteinFile = new File([proteinFileContentBlob], 'protein.pdb', {type: 'text/plain'});
            stage.loadFile(proteinFile, {defaultRepresentation: true}).then((component) => {
                obj.component = component;
            });
        },
        setView: (viewName) => {
            if (viewName === obj.views.default) {
                obj.reload();
                return;
            }
            obj.component.removeAllRepresentations();
            obj.component.addRepresentation(viewName);
        },
        render: () => {
            document.addEventListener("DOMContentLoaded", () => {
                obj.reload();
            });
        }
    }
    return obj;
}

var spinner = function () {
    const readFile = async () => {
        const contents = [];
        for (const file of document.getElementById('inputSequenceFile').files) {
            if (file) {
                const reader = new FileReader();
                reader.readAsText(file, "UTF-8");
                const fileContent = await new Promise((resolve, reject) => {
                    reader.onload = function (evt) {
                        resolve(evt.target.result);
                    }
                    reader.onerror = function (evt) {
                        reject('Error reading file')
                        document.getElementById("fileContents").innerHTML = "error reading file";
                    }
                });
                contents.push(fileContent);
            }
        }
        return contents;
    }

    $('#submitInference').on('click', function () {
        $('#submitInference').attr('disabled', true);
        $('#spinner').removeClass('invisible-spinner');

        const inferenceVal = $('#inputSequence').val();

        const socket = io();
        socket.on('connect', async function () {
            socket.emit('my event', {data: 'I\'m connected!'});
        });
    })
}

