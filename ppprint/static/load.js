// const loadingBox = document.getElementById("loading-box")
// console.log(loadingBox)
$.ajax({
    type: "GET",
    url: "/upload",
    success: function(response) {
        setTimeout(()=>{
            // loadingBox.classList.add("not-visible");
            // console.log(response)
            console.log(primary_key);
            location.href = `/load/${primary_key}`;
        }, 500)
    },
    error: function(error){
        console.log(error)
    }
})